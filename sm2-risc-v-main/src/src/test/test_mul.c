#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "fp.h"

/* ============================================================
 * 用 fp_mont_mul_ct() 作为参考实现，
 * 对拍工程中实际链接到的 fp_mont_mul()
 * ============================================================ */

static int fp_eq(const fp_t *a, const fp_t *b) {
    return (a->v[0] == b->v[0]) &&
           (a->v[1] == b->v[1]) &&
           (a->v[2] == b->v[2]) &&
           (a->v[3] == b->v[3]);
}

static int fp_cmp_local(const fp_t *a, const fp_t *b) {
    for (int i = 3; i >= 0; --i) {
        if (a->v[i] > b->v[i]) return 1;
        if (a->v[i] < b->v[i]) return -1;
    }
    return 0;
}

static void fp_print_hex(const char *name, const fp_t *a) {
    printf("%s = 0x%016llx%016llx%016llx%016llx\n",
           name,
           (unsigned long long)a->v[3],
           (unsigned long long)a->v[2],
           (unsigned long long)a->v[1],
           (unsigned long long)a->v[0]);
}

static uint64_t rand64_local(void) {
    uint64_t x = 0;
    x ^= ((uint64_t)(rand() & 0xFFFF)) <<  0;
    x ^= ((uint64_t)(rand() & 0xFFFF)) << 16;
    x ^= ((uint64_t)(rand() & 0xFFFF)) << 32;
    x ^= ((uint64_t)(rand() & 0xFFFF)) << 48;
    return x;
}

static void fp_rand_mod_p(fp_t *r) {
    do {
        r->v[0] = rand64_local();
        r->v[1] = rand64_local();
        r->v[2] = rand64_local();
        r->v[3] = rand64_local();
    } while (fp_cmp_local(r, &FP_P) >= 0);
}

static fp_t fp_make(uint64_t v0, uint64_t v1, uint64_t v2, uint64_t v3) {
    fp_t x;
    x.v[0] = v0;
    x.v[1] = v1;
    x.v[2] = v2;
    x.v[3] = v3;
    return x;
}

static inline uint64_t subb_u64_ref(uint64_t a, uint64_t b, uint64_t *borrow) {
    uint64_t c  = *borrow;
    uint64_t d  = a - b;
    uint64_t b1 = (a < b);
    uint64_t d2 = d - c;
    uint64_t b2 = (d < c);
    *borrow = b1 | b2;
    return d2;
}

static fp_t fp_p_minus_1(void) {
    fp_t x = FP_P;
    uint64_t borrow = 0;
    x.v[0] = subb_u64_ref(x.v[0], 1, &borrow);
    x.v[1] = subb_u64_ref(x.v[1], 0, &borrow);
    x.v[2] = subb_u64_ref(x.v[2], 0, &borrow);
    x.v[3] = subb_u64_ref(x.v[3], 0, &borrow);
    return x;
}

static fp_t fp_p_minus_2(void) {
    fp_t x = FP_P;
    uint64_t borrow = 0;
    x.v[0] = subb_u64_ref(x.v[0], 2, &borrow);
    x.v[1] = subb_u64_ref(x.v[1], 0, &borrow);
    x.v[2] = subb_u64_ref(x.v[2], 0, &borrow);
    x.v[3] = subb_u64_ref(x.v[3], 0, &borrow);
    return x;
}

/* 参考实现：用当前 C 版 mont_mul_ct */
static void fp_mont_mul_ref(fp_t *r, const fp_t *a, const fp_t *b) {
    fp_mont_mul_ct(r, a, b);
}

static void die_mismatch_mul(const fp_t *a, const fp_t *b,
                             const fp_t *ref, const fp_t *got) {
    printf("[FAIL] fp_mont_mul mismatch\n");
    fp_print_hex("a   ", a);
    fp_print_hex("b   ", b);
    fp_print_hex("ref ", ref);
    fp_print_hex("got ", got);
    exit(1);
}

static void test_one_pair(const fp_t *a, const fp_t *b) {
    fp_t ref_mul, got_mul;

    fp_mont_mul_ref(&ref_mul, a, b);
    fp_mont_mul(&got_mul, a, b);

    if (!fp_eq(&ref_mul, &got_mul)) {
        die_mismatch_mul(a, b, &ref_mul, &got_mul);
    }
}

static void test_properties_once(const fp_t *a, const fp_t *b, const fp_t *c) {
    fp_t x, y, z;

    /* mul(a, 0) == 0 */
    fp_mont_mul(&x, a, &FP_ZERO);
    if (!fp_eq(&x, &FP_ZERO)) {
        printf("[FAIL] property: mul(a,0) != 0\n");
        fp_print_hex("a", a);
        fp_print_hex("x", &x);
        exit(1);
    }

    /* mul(a, MONT_ONE) == a */
    fp_mont_mul(&x, a, &FP_MONT_ONE);
    if (!fp_eq(&x, a)) {
        printf("[FAIL] property: mul(a, FP_MONT_ONE) != a\n");
        fp_print_hex("a", a);
        fp_print_hex("x", &x);
        exit(1);
    }

    /* 交换律：mul(a,b) == mul(b,a) */
    fp_mont_mul(&x, a, b);
    fp_mont_mul(&y, b, a);
    if (!fp_eq(&x, &y)) {
        printf("[FAIL] property: mul(a,b) != mul(b,a)\n");
        fp_print_hex("a", a);
        fp_print_hex("b", b);
        fp_print_hex("x", &x);
        fp_print_hex("y", &y);
        exit(1);
    }

    /* 结合律：mul(mul(a,b),c) == mul(a,mul(b,c)) */
    fp_mont_mul(&x, a, b);
    fp_mont_mul(&x, &x, c);

    fp_mont_mul(&y, b, c);
    fp_mont_mul(&y, a, &y);

    if (!fp_eq(&x, &y)) {
        printf("[FAIL] property: associativity failed\n");
        fp_print_hex("a", a);
        fp_print_hex("b", b);
        fp_print_hex("c", c);
        fp_print_hex("lhs", &x);
        fp_print_hex("rhs", &y);
        exit(1);
    }

    /* 自乘与 sqr 一致：mul(a,a) == sqr(a) */
    fp_mont_mul(&x, a, a);
    fp_mont_sqr(&z, a);
    if (!fp_eq(&x, &z)) {
        printf("[FAIL] property: mul(a,a) != sqr(a)\n");
        fp_print_hex("a", a);
        fp_print_hex("mul", &x);
        fp_print_hex("sqr", &z);
        exit(1);
    }
}

static void run_edge_tests(void) {
    fp_t zero = FP_ZERO;
    fp_t one  = fp_make(1, 0, 0, 0);
    fp_t mont1 = FP_MONT_ONE;
    fp_t pm1  = fp_p_minus_1();
    fp_t pm2  = fp_p_minus_2();

    fp_t near1 = fp_make(0xffffffffffffffffULL, 0, 0, 0);
    fp_t near2 = fp_make(0xffffffffffffffffULL, 0xffffffffffffffffULL, 0, 0);
    fp_t near3 = fp_make(0xffffffffffffffffULL, 0xffffffffffffffffULL,
                         0xffffffffffffffffULL, 0);

    printf("[*] Running edge tests...\n");

    test_one_pair(&zero, &zero);
    test_one_pair(&zero, &one);
    test_one_pair(&one, &zero);
    test_one_pair(&one, &one);

    test_one_pair(&mont1, &mont1);
    test_one_pair(&mont1, &zero);
    test_one_pair(&zero, &mont1);

    test_one_pair(&pm1, &one);
    test_one_pair(&pm1, &pm1);
    test_one_pair(&pm2, &pm1);
    test_one_pair(&pm1, &pm2);

    test_one_pair(&near1, &one);
    test_one_pair(&near2, &one);
    test_one_pair(&near3, &one);

    test_one_pair(&mont1, &pm1);
    test_one_pair(&pm1, &mont1);

    /* 代数性质边界测 */
    test_properties_once(&zero, &one, &mont1);
    test_properties_once(&one, &one, &one);
    test_properties_once(&mont1, &pm1, &pm2);

    printf("[OK] Edge tests passed.\n");
}

static void run_random_tests(int rounds) {
    printf("[*] Running random tests: %d rounds...\n", rounds);

    for (int i = 0; i < rounds; ++i) {
        fp_t a, b, c;
        fp_rand_mod_p(&a);
        fp_rand_mod_p(&b);
        fp_rand_mod_p(&c);

        test_one_pair(&a, &b);
        test_properties_once(&a, &b, &c);

        if ((i + 1) % 10000 == 0) {
            printf("    passed %d\n", i + 1);
        }
    }

    printf("[OK] Random tests passed.\n");
}

int main(int argc, char **argv) {
    int rounds = 100000;

    if (argc >= 2) {
        rounds = atoi(argv[1]);
        if (rounds <= 0) {
            fprintf(stderr, "invalid rounds: %s\n", argv[1]);
            return 1;
        }
    }

    srand((unsigned)time(NULL));

#ifdef USE_RV64_ASM_MUL
    printf("[*] Testing ASM fp_mont_mul against C reference\n");
#else
    printf("[*] Testing C fp_mont_mul against C reference\n");
#endif

    run_edge_tests();
    run_random_tests(rounds);

    printf("[SUCCESS] All tests passed.\n");
    return 0;
}
