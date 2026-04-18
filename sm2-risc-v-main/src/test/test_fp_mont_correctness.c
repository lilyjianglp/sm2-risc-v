#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "fp.h"

static uint64_t g_state = 0x123456789abcdef0ULL;

static uint64_t xorshift64(void) {
    uint64_t x = g_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    g_state = x;
    return x;
}

static void rand_bytes(uint8_t *out, size_t n) {
    for (size_t i = 0; i < n; ) {
        uint64_t r = xorshift64();
        for (int j = 0; j < 8 && i < n; j++, i++) {
            out[i] = (uint8_t)(r >> (56 - 8 * j));
        }
    }
}

static void rand_fp_normal(fp_t *x) {
    uint8_t buf[32];
    rand_bytes(buf, sizeof(buf));
    fp_from_bytes(x, buf); /* 自动规约到 [0, p) */
}

static void print_fp(const char *name, const fp_t *x) {
    uint8_t buf[32];
    fp_to_bytes(buf, x);
    printf("%s = ", name);
    for (int i = 0; i < 32; i++) printf("%02x", buf[i]);
    printf("\n");
}

static int check_equal(const char *tag, const fp_t *a, const fp_t *b) {
    if (!fp_is_equal(a, b)) {
        printf("[FAIL] %s\n", tag);
        print_fp("lhs", a);
        print_fp("rhs", b);
        return 0;
    }
    return 1;
}

/* --------------------------------------------------
 * 1) round-trip:
 *    from_mont(to_mont(a)) == a
 * -------------------------------------------------- */
static int test_roundtrip_random(int rounds) {
    for (int i = 0; i < rounds; i++) {
        fp_t a, am, back;
        rand_fp_normal(&a);

        fp_to_mont(&am, &a);
        fp_from_mont(&back, &am);

        if (!check_equal("from_mont(to_mont(a)) == a", &back, &a)) {
            printf("[FAIL] roundtrip random #%d\n", i + 1);
            print_fp("a   ", &a);
            print_fp("a_m ", &am);
            print_fp("back", &back);
            return 0;
        }
    }
    return 1;
}

/* --------------------------------------------------
 * 2) to_mont preserves add:
 *    to(a+b) == mont_add(to(a), to(b))
 * -------------------------------------------------- */
static int test_add_homomorphism_random(int rounds) {
    for (int i = 0; i < rounds; i++) {
        fp_t a, b, sum;
        fp_t am, bm, lhs, rhs;

        rand_fp_normal(&a);
        rand_fp_normal(&b);

        fp_add(&sum, &a, &b);

        fp_to_mont(&lhs, &sum);

        fp_to_mont(&am, &a);
        fp_to_mont(&bm, &b);
        fp_mont_add(&rhs, &am, &bm);

        if (!check_equal("to(a+b) == mont_add(to(a),to(b))", &lhs, &rhs)) {
            printf("[FAIL] add homomorphism random #%d\n", i + 1);
            print_fp("a   ", &a);
            print_fp("b   ", &b);
            print_fp("sum ", &sum);
            print_fp("lhs ", &lhs);
            print_fp("rhs ", &rhs);
            return 0;
        }
    }
    return 1;
}

/* --------------------------------------------------
 * 3) to_mont preserves sub:
 *    to(a-b) == mont_sub(to(a), to(b))
 * -------------------------------------------------- */
static int test_sub_homomorphism_random(int rounds) {
    for (int i = 0; i < rounds; i++) {
        fp_t a, b, diff;
        fp_t am, bm, lhs, rhs;

        rand_fp_normal(&a);
        rand_fp_normal(&b);

        fp_sub(&diff, &a, &b);

        fp_to_mont(&lhs, &diff);

        fp_to_mont(&am, &a);
        fp_to_mont(&bm, &b);
        fp_mont_sub(&rhs, &am, &bm);

        if (!check_equal("to(a-b) == mont_sub(to(a),to(b))", &lhs, &rhs)) {
            printf("[FAIL] sub homomorphism random #%d\n", i + 1);
            print_fp("a    ", &a);
            print_fp("b    ", &b);
            print_fp("diff ", &diff);
            print_fp("lhs  ", &lhs);
            print_fp("rhs  ", &rhs);
            return 0;
        }
    }
    return 1;
}

/* --------------------------------------------------
 * 4) from_mont preserves add/sub:
 *    from(mont_add(x,y)) == from(x)+from(y)
 *    from(mont_sub(x,y)) == from(x)-from(y)
 * -------------------------------------------------- */
static int test_from_mont_consistency_random(int rounds) {
    for (int i = 0; i < rounds; i++) {
        fp_t a, b;
        fp_t am, bm, rm, out1, out2, ref;

        rand_fp_normal(&a);
        rand_fp_normal(&b);

        fp_to_mont(&am, &a);
        fp_to_mont(&bm, &b);

        /* add */
        fp_mont_add(&rm, &am, &bm);
        fp_from_mont(&out1, &rm);
        fp_add(&ref, &a, &b);

        if (!check_equal("from(mont_add) == add", &out1, &ref)) {
            printf("[FAIL] from_mont add consistency random #%d\n", i + 1);
            return 0;
        }

        /* sub */
        fp_mont_sub(&rm, &am, &bm);
        fp_from_mont(&out2, &rm);
        fp_sub(&ref, &a, &b);

        if (!check_equal("from(mont_sub) == sub", &out2, &ref)) {
            printf("[FAIL] from_mont sub consistency random #%d\n", i + 1);
            return 0;
        }
    }
    return 1;
}

/* --------------------------------------------------
 * 5) alias cases
 * -------------------------------------------------- */
static int test_alias_cases(int rounds) {
    for (int i = 0; i < rounds; i++) {
        fp_t a, b;
        fp_t am, bm, ref, got;

        rand_fp_normal(&a);
        rand_fp_normal(&b);

        fp_to_mont(&am, &a);
        fp_to_mont(&bm, &b);

        /* add: r == a */
        fp_mont_add(&ref, &am, &bm);
        got = am;
        fp_mont_add(&got, &got, &bm);
        if (!check_equal("alias add r==a", &ref, &got)) return 0;

        /* add: r == b */
        fp_mont_add(&ref, &am, &bm);
        got = bm;
        fp_mont_add(&got, &am, &got);
        if (!check_equal("alias add r==b", &ref, &got)) return 0;

        /* sub: r == a */
        fp_mont_sub(&ref, &am, &bm);
        got = am;
        fp_mont_sub(&got, &got, &bm);
        if (!check_equal("alias sub r==a", &ref, &got)) return 0;

        /* to_mont: r == a */
        ref = a;
        fp_to_mont(&ref, &ref);
        got = a;
        fp_to_mont(&got, &got);
        if (!check_equal("alias to_mont r==a", &ref, &got)) return 0;

        /* from_mont: r == a */
        ref = am;
        fp_from_mont(&ref, &ref);
        got = am;
        fp_from_mont(&got, &got);
        if (!check_equal("alias from_mont r==a", &ref, &got)) return 0;
    }
    return 1;
}

/* --------------------------------------------------
 * 6) special boundary cases
 * -------------------------------------------------- */
static int test_special_cases(void) {
    fp_t zero, one, pm1;
    fp_t zm, om, pm1m;
    fp_t got, ref;

    fp_set_zero(&zero);
    fp_set_one(&one);

    pm1 = FP_P;
    fp_sub(&pm1, &pm1, &one);   /* p-1 */

    fp_to_mont(&zm, &zero);
    fp_to_mont(&om, &one);
    fp_to_mont(&pm1m, &pm1);

    /* roundtrip 0 */
    fp_from_mont(&got, &zm);
    if (!check_equal("roundtrip 0", &got, &zero)) return 0;

    /* roundtrip 1 */
    fp_from_mont(&got, &om);
    if (!check_equal("roundtrip 1", &got, &one)) return 0;

    /* roundtrip p-1 */
    fp_from_mont(&got, &pm1m);
    if (!check_equal("roundtrip p-1", &got, &pm1)) return 0;

    /* (p-1)+1 = 0 */
    fp_mont_add(&got, &pm1m, &om);
    fp_from_mont(&got, &got);
    if (!check_equal("(p-1)+1 == 0", &got, &zero)) return 0;

    /* 0-1 = p-1 */
    fp_mont_sub(&got, &zm, &om);
    fp_from_mont(&got, &got);
    if (!check_equal("0-1 == p-1", &got, &pm1)) return 0;

    /* to(a+b) with a=0,b=0 */
    fp_add(&ref, &zero, &zero);
    fp_to_mont(&ref, &ref);
    fp_mont_add(&got, &zm, &zm);
    if (!check_equal("to(0+0)==mont_add", &got, &ref)) return 0;

    return 1;
}

int main(void) {
    const int N = 10000;

    printf("[*] test_fp_mont_standalone start\n");

    if (!test_special_cases()) {
        printf("[FAIL] special cases failed\n");
        return 1;
    }
    printf("[PASS] special cases\n");

    if (!test_roundtrip_random(N)) {
        printf("[FAIL] roundtrip random failed\n");
        return 1;
    }
    printf("[PASS] roundtrip random x%d\n", N);

    if (!test_add_homomorphism_random(N)) {
        printf("[FAIL] add homomorphism random failed\n");
        return 1;
    }
    printf("[PASS] add homomorphism random x%d\n", N);

    if (!test_sub_homomorphism_random(N)) {
        printf("[FAIL] sub homomorphism random failed\n");
        return 1;
    }
    printf("[PASS] sub homomorphism random x%d\n", N);

    if (!test_from_mont_consistency_random(N)) {
        printf("[FAIL] from_mont consistency random failed\n");
        return 1;
    }
    printf("[PASS] from_mont consistency random x%d\n", N);

    if (!test_alias_cases(2000)) {
        printf("[FAIL] alias cases failed\n");
        return 1;
    }
    printf("[PASS] alias cases x2000\n");

    printf("[PASS] all standalone Montgomery tests passed\n");
    return 0;
}