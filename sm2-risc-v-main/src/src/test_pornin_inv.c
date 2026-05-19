#include "fp.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef PORNIN_INV_RANDOM_TESTS
#define PORNIN_INV_RANDOM_TESTS 10000
#endif

#ifndef PORNIN_TEST_SEED
#define PORNIN_TEST_SEED 0
#endif

static void print_fp(const char *name, const fp_t *a) {
    printf("%s = 0x%016llx%016llx%016llx%016llx\n",
           name,
           (unsigned long long)a->v[3],
           (unsigned long long)a->v[2],
           (unsigned long long)a->v[1],
           (unsigned long long)a->v[0]);
}

static uint64_t rng64(void) {
    uint64_t x = 0;
    for (int i = 0; i < 4; i++) {
        x = (x << 16) ^ (uint64_t)(rand() & 0xFFFF);
    }
    return x;
}

static void fp_set_u64(fp_t *r, uint64_t x) {
    r->v[0] = x;
    r->v[1] = 0;
    r->v[2] = 0;
    r->v[3] = 0;
}

static uint64_t subb_u64_local(uint64_t a, uint64_t b, uint64_t *borrow) {
    uint64_t c = *borrow;
    uint64_t d = a - b;
    uint64_t b1 = (a < b);
    uint64_t e = d - c;
    uint64_t b2 = (d < c);
    *borrow = b1 | b2;
    return e;
}

static void fp_sub_plain_local(fp_t *r, const fp_t *a, const fp_t *b) {
    uint64_t borrow = 0;
    r->v[0] = subb_u64_local(a->v[0], b->v[0], &borrow);
    r->v[1] = subb_u64_local(a->v[1], b->v[1], &borrow);
    r->v[2] = subb_u64_local(a->v[2], b->v[2], &borrow);
    r->v[3] = subb_u64_local(a->v[3], b->v[3], &borrow);
}

#ifdef USE_PORNIN_INNER31_ASM
extern void pornin_inner_loop_31_asm(
    uint64_t xa, uint64_t xb,
    int64_t *pf0, int64_t *pg0,
    int64_t *pf1, int64_t *pg1);

static uint64_t ct_mask_u64_local(uint64_t bit01) {
    return 0ULL - (bit01 & 1ULL);
}

static uint64_t select_u64_local(uint64_t a, uint64_t b, uint64_t mask) {
    return (a & mask) | (b & ~mask);
}

static int64_t select_i64_local(int64_t a, int64_t b, uint64_t mask) {
    return (int64_t)select_u64_local((uint64_t)a, (uint64_t)b, mask);
}

static void pornin_inner_loop_31_ref(
    uint64_t xa, uint64_t xb,
    int64_t *pf0, int64_t *pg0,
    int64_t *pf1, int64_t *pg1)
{
    int64_t f0 = 1, g0 = 0;
    int64_t f1 = 0, g1 = 1;

    for (int j = 0; j < 31; ++j) {
        uint64_t odd = xa & 1ULL;
        uint64_t swap_mask = ct_mask_u64_local(odd & (uint64_t)(xa < xb));

        uint64_t old_xa = xa;
        int64_t old_f0 = f0;
        int64_t old_g0 = g0;

        xa = select_u64_local(xb, xa, swap_mask);
        xb = select_u64_local(old_xa, xb, swap_mask);
        f0 = select_i64_local(f1, f0, swap_mask);
        f1 = select_i64_local(old_f0, f1, swap_mask);
        g0 = select_i64_local(g1, g0, swap_mask);
        g1 = select_i64_local(old_g0, g1, swap_mask);

        uint64_t odd_mask = ct_mask_u64_local(odd);
        uint64_t xa_even = xa >> 1;
        uint64_t xa_odd = (xa - xb) >> 1;
        int64_t f0_odd = f0 - f1;
        int64_t g0_odd = g0 - g1;

        xa = select_u64_local(xa_odd, xa_even, odd_mask);
        f0 = select_i64_local(f0_odd, f0, odd_mask);
        g0 = select_i64_local(g0_odd, g0, odd_mask);
        f1 = f1 + f1;
        g1 = g1 + g1;
    }

    *pf0 = f0;
    *pg0 = g0;
    *pf1 = f1;
    *pg1 = g1;
}

static int test_inner31_one(uint64_t xa, uint64_t xb, const char *label) {
    int64_t cf0, cg0, cf1, cg1;
    int64_t af0, ag0, af1, ag1;

    pornin_inner_loop_31_ref(xa, xb, &cf0, &cg0, &cf1, &cg1);
    pornin_inner_loop_31_asm(xa, xb, &af0, &ag0, &af1, &ag1);

    if (cf0 != af0 || cg0 != ag0 || cf1 != af1 || cg1 != ag1) {
        printf("[FAIL] inner31 %s\n", label);
        printf("xa=0x%016llx xb=0x%016llx\n",
               (unsigned long long)xa, (unsigned long long)xb);
        printf("C   f0=%lld g0=%lld f1=%lld g1=%lld\n",
               (long long)cf0, (long long)cg0, (long long)cf1, (long long)cg1);
        printf("ASM f0=%lld g0=%lld f1=%lld g1=%lld\n",
               (long long)af0, (long long)ag0, (long long)af1, (long long)ag1);
        return 0;
    }

    return 1;
}

static int test_inner31_asm(void) {
    static const uint64_t fixed[][2] = {
        {0, 0},
        {0, 1},
        {1, 0},
        {1, 1},
        {2, 1},
        {3, 2},
        {0x7FFFFFFFULL, 0x80000000ULL},
        {0xFFFFFFFFFFFFFFFFULL, 1},
        {0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFEULL},
        {0x123456789ABCDEF0ULL, 0x0FEDCBA987654321ULL},
        {0xAAAAAAAAAAAAAAAAULL, 0x5555555555555555ULL},
    };

    for (size_t i = 0; i < sizeof(fixed) / sizeof(fixed[0]); ++i) {
        char label[64];
        snprintf(label, sizeof(label), "fixed %u", (unsigned)i);
        if (!test_inner31_one(fixed[i][0], fixed[i][1], label)) {
            return 0;
        }
    }

    for (int i = 0; i < 100000; ++i) {
        uint64_t xa = rng64();
        uint64_t xb = rng64();
        char label[64];
        snprintf(label, sizeof(label), "random %d", i);
        if (!test_inner31_one(xa, xb, label)) {
            return 0;
        }
    }

    printf("[OK] inner31 asm/c reference tests passed\n");
    return 1;
}
#endif

static void random_fp(fp_t *r);

#ifdef USE_PORNIN_UPDATE_AB_ASM
extern uint64_t pornin_lincomb_shift31_abs_asm(
    fp_t *out,
    int64_t c1, const fp_t *x,
    int64_t c2, const fp_t *y);

static uint64_t pornin_lincomb_shift31_abs_ref(
    fp_t *out,
    int64_t c1, const fp_t *x,
    int64_t c2, const fp_t *y)
{
    uint64_t z[5], mag[5];
    __int128 carry = 0;

    for (int i = 0; i < 4; ++i) {
        __int128 w = carry
            + (__int128)c1 * (__int128)x->v[i]
            + (__int128)c2 * (__int128)y->v[i];
        z[i] = (uint64_t)w;
        carry = w >> 64;
    }
    z[4] = (uint64_t)carry;

    uint64_t is_neg = (uint64_t)(carry < 0);
    uint64_t neg_mask = ct_mask_u64_local(is_neg);
    uint64_t add_one = is_neg;

    for (int i = 0; i < 5; ++i) {
        uint64_t zi = z[i] ^ neg_mask;
        uint64_t mi = zi + add_one;
        add_one = (mi < zi);
        mag[i] = mi;
    }

    out->v[0] = (mag[0] >> 31) | (mag[1] << 33);
    out->v[1] = (mag[1] >> 31) | (mag[2] << 33);
    out->v[2] = (mag[2] >> 31) | (mag[3] << 33);
    out->v[3] = (mag[3] >> 31) | (mag[4] << 33);

    return is_neg;
}

static int test_update_ab_lincomb_one(
    int64_t c1, const fp_t *x,
    int64_t c2, const fp_t *y,
    const char *label)
{
    fp_t cref, casm;
    uint64_t nref = pornin_lincomb_shift31_abs_ref(&cref, c1, x, c2, y);
    uint64_t nasm = pornin_lincomb_shift31_abs_asm(&casm, c1, x, c2, y);

    if (nref != nasm || !fp_is_equal(&cref, &casm)) {
        printf("[FAIL] update_ab lincomb %s\n", label);
        printf("c1=%lld c2=%lld ref_neg=%llu asm_neg=%llu\n",
               (long long)c1, (long long)c2,
               (unsigned long long)nref, (unsigned long long)nasm);
        print_fp("x", x);
        print_fp("y", y);
        print_fp("ref", &cref);
        print_fp("asm", &casm);
        return 0;
    }

    return 1;
}

static int64_t rand_factor31(void) {
    uint64_t u = rng64() & 0x7FFFFFFFULL;
    return (rng64() & 1ULL) ? -(int64_t)u : (int64_t)u;
}

static int test_update_ab_lincomb_asm(void) {
    fp_t x, y;
    static const int64_t factors[][2] = {
        {0, 0},
        {1, 0},
        {0, 1},
        {-1, 0},
        {0, -1},
        {1, 1},
        {-1, 1},
        {1, -1},
        {0x7FFFFFFFLL, 1},
        {-0x7FFFFFFFLL, 0x12345},
        {0x40000000LL, -0x40000000LL},
    };

    fp_set_one(&x);
    fp_set_u64(&y, 2);
    for (size_t i = 0; i < sizeof(factors) / sizeof(factors[0]); ++i) {
        char label[64];
        snprintf(label, sizeof(label), "fixed small %u", (unsigned)i);
        if (!test_update_ab_lincomb_one(factors[i][0], &x, factors[i][1], &y, label)) {
            return 0;
        }
    }

    fp_sub_plain_local(&x, &FP_P, &FP_ONE);
    y.v[0] = 0x123456789abcdef0ULL;
    y.v[1] = 0x0fedcba987654321ULL;
    y.v[2] = 0x1111222233334444ULL;
    y.v[3] = 0x5555666677778888ULL;
    for (size_t i = 0; i < sizeof(factors) / sizeof(factors[0]); ++i) {
        char label[64];
        snprintf(label, sizeof(label), "fixed wide %u", (unsigned)i);
        if (!test_update_ab_lincomb_one(factors[i][0], &x, factors[i][1], &y, label)) {
            return 0;
        }
    }

    for (int i = 0; i < 100000; ++i) {
        random_fp(&x);
        random_fp(&y);
        char label[64];
        snprintf(label, sizeof(label), "random %d", i);
        if (!test_update_ab_lincomb_one(rand_factor31(), &x, rand_factor31(), &y, label)) {
            return 0;
        }
    }

    printf("[OK] update_ab lincomb asm/c reference tests passed\n");
    return 1;
}
#endif

static void random_fp(fp_t *r) {
    r->v[0] = rng64();
    r->v[1] = rng64();
    r->v[2] = rng64();
    r->v[3] = rng64();

    /*
     * SM2 p 很接近 2^256，所以随机 256-bit 数减一次 p 即可落到 [0,p)。
     */
    if (fp_cmp(r, &FP_P) >= 0) {
        fp_sub_plain_local(r, r, &FP_P);
    }

    if (fp_is_zero(r)) {
        fp_set_one(r);
    }
}

/*
 * 用 Montgomery 乘法临时验证普通域乘法：
 *   r = a*b mod p
 */
static void fp_mul_plain_via_mont(fp_t *r, const fp_t *a, const fp_t *b) {
    fp_t am, bm, rm;

    fp_to_mont(&am, a);
    fp_to_mont(&bm, b);
    fp_mont_mul(&rm, &am, &bm);
    fp_from_mont(r, &rm);
}

static int test_one(const fp_t *a, const char *label) {
    fp_t inv, prod;

    fp_pornin_inv(&inv, a);
    fp_mul_plain_via_mont(&prod, a, &inv);

    if (!fp_is_equal(&prod, &FP_ONE)) {
        printf("[FAIL] %s\n", label);
        print_fp("a", a);
        print_fp("inv", &inv);
        print_fp("a*inv", &prod);
        return 0;
    }

    return 1;
}

static int test_one_mont(const fp_t *a, const char *label) {
    fp_t a_bar, inv_bar, prod_bar, prod;

    fp_to_mont(&a_bar, a);
    fp_mont_inv(&inv_bar, &a_bar);
    fp_mont_mul(&prod_bar, &a_bar, &inv_bar);
    fp_from_mont(&prod, &prod_bar);

    if (!fp_is_equal(&prod, &FP_ONE)) {
        printf("[FAIL] mont wrapper %s\n", label);
        print_fp("a", a);
        print_fp("a_bar", &a_bar);
        print_fp("inv_bar", &inv_bar);
        print_fp("a*inv", &prod);
        return 0;
    }

    return 1;
}

static int test_both(const fp_t *a, const char *label) {
    return test_one(a, label) && test_one_mont(a, label);
}

int main(void) {
    unsigned seed = PORNIN_TEST_SEED ? (unsigned)PORNIN_TEST_SEED : (unsigned)time(NULL);
    srand(seed);

#ifdef PORNIN_PROFILE
    pornin_profile_reset();
#endif

#ifdef USE_PORNIN_INNER31_ASM
    if (!test_inner31_asm()) return 1;
#endif
#ifdef USE_PORNIN_UPDATE_AB_ASM
    if (!test_update_ab_lincomb_asm()) return 1;
#endif

    fp_t a;

    /*
     * 固定测试 1：a = 1
     */
    fp_set_one(&a);
    if (!test_both(&a, "a = 1")) return 1;

    /*
     * 固定测试 2：a = 2
     */
    fp_set_u64(&a, 2);
    if (!test_both(&a, "a = 2")) return 1;

    /*
     * 固定测试 3：a = 3
     */
    fp_set_u64(&a, 3);
    if (!test_both(&a, "a = 3")) return 1;

    /*
     * 固定测试 4：a = p - 1
     */
    fp_sub_plain_local(&a, &FP_P, &FP_ONE);
    if (!test_both(&a, "a = p-1")) return 1;

    /*
     * 固定测试 5：普通随机样例
     */
    a.v[0] = 0x123456789abcdef0ULL;
    a.v[1] = 0x0fedcba987654321ULL;
    a.v[2] = 0x1111222233334444ULL;
    a.v[3] = 0x5555666677778888ULL;
    if (fp_cmp(&a, &FP_P) >= 0) {
        fp_sub_plain_local(&a, &a, &FP_P);
    }
    if (!test_both(&a, "fixed random-looking")) return 1;

    printf("[OK] fixed tests passed\n");

    /*
     * 随机测试
     */
    printf("[INFO] random seed=%u random tests=%d\n", seed, PORNIN_INV_RANDOM_TESTS);
    for (int i = 0; i < PORNIN_INV_RANDOM_TESTS; i++) {
        random_fp(&a);

        char label[64];
        snprintf(label, sizeof(label), "random %d", i);

        if (!test_both(&a, label)) {
            return 1;
        }
    }

    printf("[OK] random tests passed\n");

    /*
     * 测试 a = 0，按你的实现应该返回 0。
     */
    a.v[0] = 0;
    a.v[1] = 0;
    a.v[2] = 0;
    a.v[3] = 0;

    fp_t inv0;
    fp_pornin_inv(&inv0, &a);

    if (!fp_is_zero(&inv0)) {
        printf("[FAIL] zero input should return zero\n");
        print_fp("inv0", &inv0);
        return 1;
    }

    fp_t zero_bar, inv0_bar;
    fp_to_mont(&zero_bar, &a);
    fp_mont_inv(&inv0_bar, &zero_bar);

    if (!fp_is_zero(&inv0_bar)) {
        printf("[FAIL] mont zero input should return zero\n");
        print_fp("inv0_bar", &inv0_bar);
        return 1;
    }

    printf("[OK] zero test passed\n");
#ifdef PORNIN_PROFILE
    pornin_profile_print();
#endif
    printf("all tests passed\n");

    return 0;
}
