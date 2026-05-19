/*
 * fp_pornin_inv_sm2_16x31_14.c
 *
 * SM2 p = FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF
 * Pornin optimized binary GCD inversion, 16*31 + 14 = 510 iterations.
 *
 * This file exports ordinary-domain inversion:
 *      void fp_pornin_inv(fp_t *r, const fp_t *a)
 *
 * Semantics:
 *      input : a in [0, p)
 *      output: a^{-1} mod p, or 0 if a == 0
 *
 * Notes:
 *   - The inner/tail divstep loops use mask-based selects instead of
 *     data-dependent swaps/branches.
 *   - It reuses existing Montgomery multiplication only for the final multiply
 *     by C510 = 2^{-510} mod p.
 */

#include "fp.h"
#include <stdint.h>
#ifdef PORNIN_PROFILE
#include <stdio.h>
#include <time.h>
#endif

#ifdef PORNIN_PERF_NOINLINE
#define PORNIN_NOINLINE __attribute__((noinline))
#else
#define PORNIN_NOINLINE
#endif

/* C510 in Montgomery form: C510*R mod p. */
static const fp_t FP_INV2_POW_510_BAR = {
    .v = {
        0xFFFFFFE400000013ULL,
        0xFFFFFFF700000018ULL,
        0xFFFFFFF00000000BULL,
        0xFFFFFFEF00000017ULL
    }
};

/* R^2 mod p, used to return a Montgomery-domain inverse directly. */
static const fp_t FP_R2 = {
    .v = {
        0x0000000200000003ULL,
        0x00000002FFFFFFFFULL,
        0x0000000100000001ULL,
        0x0000000400000002ULL
    }
};

#ifdef USE_PORNIN_FULL_INV_ASM
const fp_t PORNIN_FP_INV2_POW_510_BAR = {
    .v = {
        0xFFFFFFE400000013ULL,
        0xFFFFFFF700000018ULL,
        0xFFFFFFF00000000BULL,
        0xFFFFFFEF00000017ULL
    }
};

const fp_t PORNIN_FP_R2 = {
    .v = {
        0x0000000200000003ULL,
        0x00000002FFFFFFFFULL,
        0x0000000100000001ULL,
        0x0000000400000002ULL
    }
};
#endif

#ifdef PORNIN_PROFILE
typedef struct {
    uint64_t calls;
    uint64_t zero_calls;
    uint64_t cycle_probe;
    uint64_t rounds;
    uint64_t extract;
    uint64_t inner31;
    uint64_t update_ab;
    uint64_t update_uv;
    uint64_t tail14;
    uint64_t tail_uv;
    uint64_t final_mul;
} pornin_profile_t;

static pornin_profile_t g_pornin_profile;

static __attribute__((noinline)) uint64_t pornin_rdcycle(void)
{
#if defined(PORNIN_USE_RDTIME) && defined(__riscv)
    uint64_t x;
    __asm__ volatile ("" ::: "memory");
    __asm__ volatile ("csrr %0, time" : "=r"(x) :: "memory");
    __asm__ volatile ("" ::: "memory");
    return x;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
#endif
}

void pornin_profile_reset(void)
{
    g_pornin_profile = (pornin_profile_t){0};
    uint64_t t0 = pornin_rdcycle();
    uint64_t t1 = pornin_rdcycle();
    g_pornin_profile.cycle_probe = t1 - t0;
}

static void pornin_profile_print_line(const char *name, uint64_t cycles, uint64_t denom)
{
    if (denom == 0) {
        printf("%-18s : n/a\n", name);
        return;
    }
    printf("%-18s : total=%llu ns avg=%llu ns\n",
           name,
           (unsigned long long)cycles,
           (unsigned long long)(cycles / denom));
}

void pornin_profile_print(void)
{
    uint64_t calls = g_pornin_profile.calls;
    uint64_t rounds = g_pornin_profile.rounds;

    printf("\n=== Pornin inversion profile ===\n");
    printf("calls              : %llu\n", (unsigned long long)calls);
    printf("zero calls         : %llu\n", (unsigned long long)g_pornin_profile.zero_calls);
    printf("timer probe        : %llu\n", (unsigned long long)g_pornin_profile.cycle_probe);
    printf("rounds             : %llu\n", (unsigned long long)rounds);
    pornin_profile_print_line("extract", g_pornin_profile.extract, rounds);
    pornin_profile_print_line("inner31", g_pornin_profile.inner31, rounds);
    pornin_profile_print_line("update_ab", g_pornin_profile.update_ab, rounds);
    pornin_profile_print_line("update_uv", g_pornin_profile.update_uv, rounds);
    pornin_profile_print_line("tail14", g_pornin_profile.tail14, calls);
    pornin_profile_print_line("tail_uv", g_pornin_profile.tail_uv, calls);
    pornin_profile_print_line("final_mul", g_pornin_profile.final_mul, calls);
    printf("================================\n");
}
#endif

static inline uint64_t pornin_ct_mask_u64(uint64_t bit01)
{
    return 0ULL - (bit01 & 1ULL);
}

static inline uint64_t pornin_select_u64(uint64_t a, uint64_t b, uint64_t mask)
{
    return (a & mask) | (b & ~mask);
}

static inline int64_t pornin_select_i64(int64_t a, int64_t b, uint64_t mask)
{
    uint64_t ua = (uint64_t)a;
    uint64_t ub = (uint64_t)b;
    return (int64_t)pornin_select_u64(ua, ub, mask);
}

static int pornin_u256_bitlen(const fp_t *a)
{
    uint64_t seen = 0;
    int bitlen = 0;

    for (int i = 3; i >= 0; --i) {
        uint64_t x = a->v[i];
        uint64_t nz = (uint64_t)(x != 0);
        uint64_t take = nz & (seen ^ 1ULL);
        uint64_t take_mask = pornin_ct_mask_u64(take);
        int limb_bits;

#if defined(__GNUC__) || defined(__clang__)
        limb_bits = 64 - __builtin_clzll(x | (nz ^ 1ULL));
#else
        uint64_t y = x;
        limb_bits = 0;
        for (int j = 0; j < 64; ++j) {
            limb_bits += (int)(y != 0);
            y >>= 1;
        }
#endif

        int candidate = i * 64 + limb_bits;
        bitlen = (int)(((uint64_t)candidate & take_mask) |
                       ((uint64_t)bitlen & ~take_mask));
        seen |= nz;
    }

    return bitlen;
}

static uint64_t pornin_u256_shr_low64(const fp_t *x, unsigned shift)
{
    unsigned limb = shift >> 6;
    unsigned off  = shift & 63;

    uint64_t lo =
        (x->v[0] & pornin_ct_mask_u64((uint64_t)(limb == 0))) |
        (x->v[1] & pornin_ct_mask_u64((uint64_t)(limb == 1))) |
        (x->v[2] & pornin_ct_mask_u64((uint64_t)(limb == 2))) |
        (x->v[3] & pornin_ct_mask_u64((uint64_t)(limb == 3)));

    uint64_t hi =
        (x->v[1] & pornin_ct_mask_u64((uint64_t)(limb == 0))) |
        (x->v[2] & pornin_ct_mask_u64((uint64_t)(limb == 1))) |
        (x->v[3] & pornin_ct_mask_u64((uint64_t)(limb == 2)));

    uint64_t hi_part = hi << ((64 - off) & 63);
    hi_part &= pornin_ct_mask_u64((uint64_t)(off != 0));

    return (lo >> off) | hi_part;
}

/*
 * Build 64-bit approximations for k = 32:
 *   x_hat = low31(x) | (high33(x) << 31)
 * with high33(x) = floor(x / 2^(n - 33)).
 */
PORNIN_NOINLINE static void pornin_extract_approx64(
    const fp_t *a, const fp_t *b,
    uint64_t *a_hat, uint64_t *b_hat)
{
    int na = pornin_u256_bitlen(a);
    int nb = pornin_u256_bitlen(b);
    int n  = (na > nb) ? na : nb;

    uint64_t n_lt_64_mask = pornin_ct_mask_u64((uint64_t)(n < 64));
    n = (int)(((uint64_t)64 & n_lt_64_mask) | ((uint64_t)n & ~n_lt_64_mask));

    /* Correct shift is n - 33, not n - 64. */
    unsigned shift = (unsigned)(n - 33);

    uint64_t hi_a = pornin_u256_shr_low64(a, shift) & 0x1FFFFFFFFULL;
    uint64_t hi_b = pornin_u256_shr_low64(b, shift) & 0x1FFFFFFFFULL;
    uint64_t lo_a = a->v[0] & 0x7FFFFFFFULL;
    uint64_t lo_b = b->v[0] & 0x7FFFFFFFULL;

    *a_hat = lo_a | (hi_a << 31);
    *b_hat = lo_b | (hi_b << 31);
}

/* Mask-based inner loop: 31 iterations.
 *
 * Build option:
 *   - default: use the portable C implementation below
 *   - with -DUSE_PORNIN_INNER31_ASM:
 *       call asm/pornin_inner31.S:
 *           pornin_inner_loop_31_asm(xa, xb, pf0, pg0, pf1, pg1)
 */
#ifdef USE_PORNIN_INNER31_ASM
extern void pornin_inner_loop_31_asm(
    uint64_t xa, uint64_t xb,
    int64_t *pf0, int64_t *pg0,
    int64_t *pf1, int64_t *pg1);

static void pornin_inner_loop_31(
    uint64_t xa, uint64_t xb,
    int64_t *pf0, int64_t *pg0,
    int64_t *pf1, int64_t *pg1)
{
    pornin_inner_loop_31_asm(xa, xb, pf0, pg0, pf1, pg1);
}
#else
static void pornin_inner_loop_31(
    uint64_t xa, uint64_t xb,
    int64_t *pf0, int64_t *pg0,
    int64_t *pf1, int64_t *pg1)
{
    int64_t f0 = 1, g0 = 0;
    int64_t f1 = 0, g1 = 1;

    for (int j = 0; j < 31; ++j) {
        uint64_t odd = xa & 1ULL;
        uint64_t swap_mask = pornin_ct_mask_u64(odd & (uint64_t)(xa < xb));

        uint64_t txa = xa;
        int64_t tf0 = f0;
        int64_t tg0 = g0;

        xa = pornin_select_u64(xb, xa, swap_mask);
        xb = pornin_select_u64(txa, xb, swap_mask);
        f0 = pornin_select_i64(f1, f0, swap_mask);
        f1 = pornin_select_i64(tf0, f1, swap_mask);
        g0 = pornin_select_i64(g1, g0, swap_mask);
        g1 = pornin_select_i64(tg0, g1, swap_mask);

        uint64_t odd_mask = pornin_ct_mask_u64(odd);
        uint64_t xa_even = xa >> 1;
        uint64_t xa_odd = (xa - xb) >> 1;
        int64_t f0_odd = f0 - f1;
        int64_t g0_odd = g0 - g1;

        xa = pornin_select_u64(xa_odd, xa_even, odd_mask);
        f0 = pornin_select_i64(f0_odd, f0, odd_mask);
        g0 = pornin_select_i64(g0_odd, g0, odd_mask);
        f1 = f1 + f1;
        g1 = g1 + g1;
    }

    *pf0 = f0;
    *pg0 = g0;
    *pf1 = f1;
    *pg1 = g1;
}
#endif

/* ============================================================
 *  320-bit helpers for abs((c1*x + c2*y) >> 31)
 * ============================================================ */

/*
 * out = abs((c1*x + c2*y) >> 31)
 * Return 1 if the original signed value was negative, else 0.
 *
 * This follows the s256_lin_div31_abs idea: compute one signed 320-bit
 * linear combination, then take the absolute value in two's-complement
 * form before the exact shift.
 */
PORNIN_NOINLINE static uint64_t pornin_lincomb_shift31_abs_c(
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
    uint64_t neg_mask = pornin_ct_mask_u64(is_neg);
    uint64_t add_one = is_neg;

    for (int i = 0; i < 5; ++i) {
        uint64_t zi = z[i] ^ neg_mask;
        uint64_t mi = zi + add_one;
        add_one = (mi < zi);
        mag[i] = mi;
    }

    /* Exact right shift by 31. */
    out->v[0] = (mag[0] >> 31) | (mag[1] << 33);
    out->v[1] = (mag[1] >> 31) | (mag[2] << 33);
    out->v[2] = (mag[2] >> 31) | (mag[3] << 33);
    out->v[3] = (mag[3] >> 31) | (mag[4] << 33);

    return is_neg;
}

#define pornin_lincomb_shift31_abs pornin_lincomb_shift31_abs_c

PORNIN_NOINLINE static void pornin_update_ab_c(
    fp_t *new_a, fp_t *new_b,
    const fp_t *a, const fp_t *b,
    int64_t *f0, int64_t *g0,
    int64_t *f1, int64_t *g1)
{
    uint64_t neg_a = pornin_lincomb_shift31_abs(new_a, *f0, a, *g0, b);
    uint64_t neg_a_mask = pornin_ct_mask_u64(neg_a);
    *f0 = pornin_select_i64(-*f0, *f0, neg_a_mask);
    *g0 = pornin_select_i64(-*g0, *g0, neg_a_mask);

    uint64_t neg_b = pornin_lincomb_shift31_abs(new_b, *f1, a, *g1, b);
    uint64_t neg_b_mask = pornin_ct_mask_u64(neg_b);
    *f1 = pornin_select_i64(-*f1, *f1, neg_b_mask);
    *g1 = pornin_select_i64(-*g1, *g1, neg_b_mask);
}

#ifdef USE_PORNIN_UPDATE_AB_ASM
extern void pornin_update_ab_asm(
    fp_t *new_a, fp_t *new_b,
    const fp_t *a, const fp_t *b,
    int64_t *f0, int64_t *g0,
    int64_t *f1, int64_t *g1);
#define pornin_update_ab pornin_update_ab_asm
#else
#define pornin_update_ab pornin_update_ab_c
#endif

/* ============================================================
 *  Fast small signed linear combination modulo SM2 p.
 *
 *      r = c1*x + c2*y mod p
 *
 *  This replaces the old 32-step double-and-add version.
 *
 *  Method:
 *    1. If coefficient is negative, replace x by -x mod p.
 *    2. Compute abs(c1)*x + abs(c2)*y into 5 limbs.
 *    3. Reduce with SM2 folding:
 *
 *       2^256 = 2^224 + 2^96 - 2^64 + 1 mod p
 * ============================================================ */

static inline uint64_t pornin_subb_local(uint64_t a, uint64_t b, uint64_t *borrow)
{
    uint64_t c = *borrow;
    uint64_t d = a - b;
    uint64_t b1 = (a < b);
    uint64_t d2 = d - c;
    uint64_t b2 = (d < c);
    *borrow = b1 | b2;
    return d2;
}

static int pornin_fp_is_zero_local(const fp_t *a)
{
    return (a->v[0] | a->v[1] | a->v[2] | a->v[3]) == 0;
}

static void pornin_fp_neg_mod(fp_t *r, const fp_t *a)
{
    fp_t t;
    uint64_t borrow = 0;

    t.v[0] = pornin_subb_local(FP_P.v[0], a->v[0], &borrow);
    t.v[1] = pornin_subb_local(FP_P.v[1], a->v[1], &borrow);
    t.v[2] = pornin_subb_local(FP_P.v[2], a->v[2], &borrow);
    t.v[3] = pornin_subb_local(FP_P.v[3], a->v[3], &borrow);

    /*
     * If a == 0, return 0 instead of p.
     */
    uint64_t zero_mask = pornin_ct_mask_u64((uint64_t)pornin_fp_is_zero_local(a));
    uint64_t nonzero_mask = ~zero_mask;

    r->v[0] = t.v[0] & nonzero_mask;
    r->v[1] = t.v[1] & nonzero_mask;
    r->v[2] = t.v[2] & nonzero_mask;
    r->v[3] = t.v[3] & nonzero_mask;
}

static void pornin_fp_cond_sub_p_local(fp_t *r)
{
    fp_t t;
    uint64_t borrow = 0;

    t.v[0] = pornin_subb_local(r->v[0], FP_P.v[0], &borrow);
    t.v[1] = pornin_subb_local(r->v[1], FP_P.v[1], &borrow);
    t.v[2] = pornin_subb_local(r->v[2], FP_P.v[2], &borrow);
    t.v[3] = pornin_subb_local(r->v[3], FP_P.v[3], &borrow);

    uint64_t m = pornin_ct_mask_u64(borrow ^ 1ULL);

    r->v[0] = (t.v[0] & m) | (r->v[0] & ~m);
    r->v[1] = (t.v[1] & m) | (r->v[1] & ~m);
    r->v[2] = (t.v[2] & m) | (r->v[2] & ~m);
    r->v[3] = (t.v[3] & m) | (r->v[3] & ~m);
}

/*
 * t += c * x
 *
 * t is 5 limbs.
 * c is small non-negative integer.
 * x is 256-bit.
 */
static void pornin_u320_addmul_small(uint64_t t[5], uint64_t c, const fp_t *x)
{
    __uint128_t carry = 0;

    for (int i = 0; i < 4; i++) {
        __uint128_t z =
            (__uint128_t)c * (__uint128_t)x->v[i]
            + (__uint128_t)t[i]
            + carry;

        t[i] = (uint64_t)z;
        carry = z >> 64;
    }

    t[4] += (uint64_t)carry;
}

/*
 * Fold a 5-limb non-negative integer modulo SM2 p.
 *
 * Input:
 *   X = t0 + t1*2^64 + t2*2^128 + t3*2^192 + h*2^256
 *
 * SM2:
 *   2^256 = 2^224 + 2^96 - 2^64 + 1 mod p
 */
PORNIN_NOINLINE static void pornin_sm2_fold_u320(fp_t *r, const uint64_t t[5])
{
    uint64_t x0 = t[0];
    uint64_t x1 = t[1];
    uint64_t x2 = t[2];
    uint64_t x3 = t[3];
    uint64_t h  = t[4];

    /*
     * The small coefficients produced by the 31-step matrix keep h modest.
     * Two fixed folds are enough to absorb the carry generated by folding
     * h*2^256 through the SM2 pseudo-Mersenne relation.
     */
    for (int round = 0; round < 2; round++) {
        __uint128_t a0 = x0;
        __uint128_t a1 = x1;
        __uint128_t a2 = x2;
        __uint128_t a3 = x3;
        __uint128_t a4 = 0;

        /*
         * + h
         */
        a0 += h;

        /*
         * + h*2^96 - h*2^64
         *
         * h*2^96:
         *   limb1 += h << 32
         *   limb2 += h >> 32
         *
         * -h*2^64:
         *   limb1 -= h
         *
         * Since h << 32 >= h for h >= 0, this expression is non-negative.
         */
        a1 += (((__uint128_t)h) << 32) - (__uint128_t)h;
        a2 += h >> 32;

        /*
         * + h*2^224
         *
         * limb3 += h << 32
         * high  += h >> 32
         */
        a3 += ((__uint128_t)h) << 32;
        a4 += h >> 32;

        x0 = (uint64_t)a0;
        a1 += a0 >> 64;

        x1 = (uint64_t)a1;
        a2 += a1 >> 64;

        x2 = (uint64_t)a2;
        a3 += a2 >> 64;

        x3 = (uint64_t)a3;
        a4 += a3 >> 64;

        h = (uint64_t)a4;
    }

    r->v[0] = x0;
    r->v[1] = x1;
    r->v[2] = x2;
    r->v[3] = x3;

    /*
     * Bring close representative into [0,p).
     */
    pornin_fp_cond_sub_p_local(r);
    pornin_fp_cond_sub_p_local(r);
    pornin_fp_cond_sub_p_local(r);
}

#ifdef USE_PORNIN_LINCOMB_ASM
void pornin_sm2_fold_u320_public(fp_t *r, const uint64_t t[5])
{
    pornin_sm2_fold_u320(r, t);
}

extern void pornin_fp_lincomb_small_mod_precomp_asm(
    fp_t *r,
    int64_t c1, const fp_t *x, const fp_t *neg_x,
    int64_t c2, const fp_t *y, const fp_t *neg_y);

#define pornin_fp_lincomb_small_mod_precomp pornin_fp_lincomb_small_mod_precomp_asm
#else
PORNIN_NOINLINE static void pornin_fp_lincomb_small_mod_precomp(
    fp_t *r,
    int64_t c1, const fp_t *x, const fp_t *neg_x,
    int64_t c2, const fp_t *y, const fp_t *neg_y)
{
    uint64_t t[5] = {0, 0, 0, 0, 0};

    uint64_t uc1 = (uint64_t)c1;
    uint64_t neg1 = uc1 >> 63;
    uint64_t m1 = 0ULL - neg1;
    uint64_t ac1 = (uc1 ^ m1) - m1;

    uint64_t uc2 = (uint64_t)c2;
    uint64_t neg2 = uc2 >> 63;
    uint64_t m2 = 0ULL - neg2;
    uint64_t ac2 = (uc2 ^ m2) - m2;

    fp_t sx;
    fp_t sy;

    uint64_t neg1_mask = pornin_ct_mask_u64(neg1);
    sx.v[0] = pornin_select_u64(neg_x->v[0], x->v[0], neg1_mask);
    sx.v[1] = pornin_select_u64(neg_x->v[1], x->v[1], neg1_mask);
    sx.v[2] = pornin_select_u64(neg_x->v[2], x->v[2], neg1_mask);
    sx.v[3] = pornin_select_u64(neg_x->v[3], x->v[3], neg1_mask);

    uint64_t neg2_mask = pornin_ct_mask_u64(neg2);
    sy.v[0] = pornin_select_u64(neg_y->v[0], y->v[0], neg2_mask);
    sy.v[1] = pornin_select_u64(neg_y->v[1], y->v[1], neg2_mask);
    sy.v[2] = pornin_select_u64(neg_y->v[2], y->v[2], neg2_mask);
    sy.v[3] = pornin_select_u64(neg_y->v[3], y->v[3], neg2_mask);

    pornin_u320_addmul_small(t, ac1, &sx);
    pornin_u320_addmul_small(t, ac2, &sy);

    pornin_sm2_fold_u320(r, t);
}
#endif

/*
 * r = c1*x + c2*y mod p
 */
static void pornin_fp_lincomb_small_mod(
    fp_t *r,
    int64_t c1, const fp_t *x,
    int64_t c2, const fp_t *y)
{
    fp_t neg_x, neg_y;
    pornin_fp_neg_mod(&neg_x, x);
    pornin_fp_neg_mod(&neg_y, y);
    pornin_fp_lincomb_small_mod_precomp(r, c1, x, &neg_x, c2, y, &neg_y);
}

#ifdef USE_PORNIN_UPDATE_UV_ASM
extern void pornin_update_uv_asm(
    fp_t *new_u, fp_t *new_v,
    const fp_t *u, const fp_t *v,
    int64_t f0, int64_t g0,
    int64_t f1, int64_t g1);

#define pornin_update_uv pornin_update_uv_asm
#else
PORNIN_NOINLINE static void pornin_update_uv(
    fp_t *new_u, fp_t *new_v,
    const fp_t *u, const fp_t *v,
    int64_t f0, int64_t g0,
    int64_t f1, int64_t g1)
{
    fp_t ru, rv, neg_u, neg_v;

    pornin_fp_neg_mod(&neg_u, u);
    pornin_fp_neg_mod(&neg_v, v);
    pornin_fp_lincomb_small_mod_precomp(&ru, f0, u, &neg_u, g0, v, &neg_v);
    pornin_fp_lincomb_small_mod_precomp(&rv, f1, u, &neg_u, g1, v, &neg_v);
    *new_u = ru;
    *new_v = rv;
}
#endif

/* Tail: 14 exact mask-based iterations on small real a,b; return second row F1,G1. */
PORNIN_NOINLINE static void pornin_tail_loop_14(const fp_t *a, const fp_t *b, int64_t *pF1, int64_t *pG1)
{
    uint64_t alpha = a->v[0];
    uint64_t beta  = b->v[0];

    int64_t F0 = 1, G0 = 0;
    int64_t F1 = 0, G1 = 1;

    for (int j = 0; j < 14; ++j) {
        uint64_t odd = alpha & 1ULL;
        uint64_t swap_mask = pornin_ct_mask_u64(odd & (uint64_t)(alpha < beta));

        uint64_t old_alpha = alpha;
        int64_t old_F0 = F0;
        int64_t old_G0 = G0;

        alpha = pornin_select_u64(beta, alpha, swap_mask);
        beta = pornin_select_u64(old_alpha, beta, swap_mask);
        F0 = pornin_select_i64(F1, F0, swap_mask);
        F1 = pornin_select_i64(old_F0, F1, swap_mask);
        G0 = pornin_select_i64(G1, G0, swap_mask);
        G1 = pornin_select_i64(old_G0, G1, swap_mask);

        uint64_t odd_mask = pornin_ct_mask_u64(odd);
        uint64_t alpha_even = alpha >> 1;
        uint64_t alpha_odd = (alpha - beta) >> 1;
        int64_t F0_odd = F0 - F1;
        int64_t G0_odd = G0 - G1;

        alpha = pornin_select_u64(alpha_odd, alpha_even, odd_mask);
        F0 = pornin_select_i64(F0_odd, F0, odd_mask);
        G0 = pornin_select_i64(G0_odd, G0, odd_mask);
        F1 = F1 + F1;
        G1 = G1 + G1;
    }

    *pF1 = F1;
    *pG1 = G1;
}

static void pornin_inv_scaled(fp_t *r, const fp_t *y, const fp_t *u0)
{
    if (fp_is_zero(y)) {
#ifdef PORNIN_PROFILE
        g_pornin_profile.zero_calls++;
#endif
        fp_set_zero(r);
        return;
    }

#ifdef PORNIN_PROFILE
    g_pornin_profile.calls++;
#endif

    fp_t a = *y;
    fp_t b = FP_P;
    fp_t u = *u0;
    fp_t v = FP_ZERO;

    for (int i = 0; i < 16; ++i) {
        uint64_t a_hat, b_hat;
#ifdef PORNIN_PROFILE
        uint64_t t0 = pornin_rdcycle();
#endif
        pornin_extract_approx64(&a, &b, &a_hat, &b_hat);
#ifdef PORNIN_PROFILE
        uint64_t t1 = pornin_rdcycle();
        g_pornin_profile.extract += t1 - t0;
#endif

        int64_t f0, g0, f1, g1;
        pornin_inner_loop_31(a_hat, b_hat, &f0, &g0, &f1, &g1);
#ifdef PORNIN_PROFILE
        uint64_t t2 = pornin_rdcycle();
        g_pornin_profile.inner31 += t2 - t1;
#endif

        fp_t new_a, new_b;
        pornin_update_ab(&new_a, &new_b, &a, &b, &f0, &g0, &f1, &g1);
#ifdef PORNIN_PROFILE
        uint64_t t3 = pornin_rdcycle();
        g_pornin_profile.update_ab += t3 - t2;
#endif

        fp_t new_u, new_v;
        pornin_update_uv(&new_u, &new_v, &u, &v, f0, g0, f1, g1);
#ifdef PORNIN_PROFILE
        uint64_t t4 = pornin_rdcycle();
        g_pornin_profile.update_uv += t4 - t3;
        g_pornin_profile.rounds++;
#endif

        a = new_a;
        b = new_b;
        u = new_u;
        v = new_v;
    }

    int64_t F1, G1;
#ifdef PORNIN_PROFILE
    uint64_t t0 = pornin_rdcycle();
#endif
    pornin_tail_loop_14(&a, &b, &F1, &G1);
#ifdef PORNIN_PROFILE
    uint64_t t1 = pornin_rdcycle();
    g_pornin_profile.tail14 += t1 - t0;
#endif

    /* Official-style tail: v = F1*u + G1*v mod p. */
    fp_t new_v;
    pornin_fp_lincomb_small_mod(&new_v, F1, &u, G1, &v);
#ifdef PORNIN_PROFILE
    uint64_t t2 = pornin_rdcycle();
    g_pornin_profile.tail_uv += t2 - t1;
#endif
    v = new_v;

    /*
     * v = v * C510 mod p.
     *
     * Montgomery multiplication computes x*y*R^-1.  Since C510_BAR is
     * C510*R mod p, MontMul(v, C510_BAR) returns v*C510 in the ordinary
     * field domain even though v itself is not in Montgomery form.
     */
    fp_mont_mul(r, &v, &FP_INV2_POW_510_BAR);
#ifdef PORNIN_PROFILE
    uint64_t t3 = pornin_rdcycle();
    g_pornin_profile.final_mul += t3 - t2;
#endif
}

static void pornin_inv_plain(fp_t *r, const fp_t *y)
{
    pornin_inv_scaled(r, y, &FP_ONE);
}

#ifdef USE_PORNIN_FULL_INV_ASM
void pornin_extract_approx64_public(
    const fp_t *a, const fp_t *b,
    uint64_t *a_hat, uint64_t *b_hat)
{
    pornin_extract_approx64(a, b, a_hat, b_hat);
}

void pornin_tail_loop_14_public(
    const fp_t *a, const fp_t *b,
    int64_t *pF1, int64_t *pG1)
{
    pornin_tail_loop_14(a, b, pF1, pG1);
}

void pornin_fp_lincomb_small_mod_public(
    fp_t *r,
    int64_t c1, const fp_t *x,
    int64_t c2, const fp_t *y)
{
    pornin_fp_lincomb_small_mod(r, c1, x, c2, y);
}

extern void fp_mont_inv_pornin_asm(fp_t *r, const fp_t *a_bar);
#endif

void fp_pornin_inv(fp_t *r, const fp_t *a)
{
    pornin_inv_plain(r, a);
}

#ifdef USE_PORNIN_INV
void fp_mont_inv(fp_t *r, const fp_t *a_bar)
{
#ifdef USE_PORNIN_FULL_INV_ASM
    fp_mont_inv_pornin_asm(r, a_bar);
#else
    /*
     * a_bar = a*R mod p.  pornin_inv_scaled() returns u0*y^(-1) mod p.
     * With u0 = R^2 and y = a_bar:
     *     R^2 * (aR)^(-1) = a^(-1)R mod p
     * which is the Montgomery-domain inverse expected by callers.
     */
    pornin_inv_scaled(r, a_bar, &FP_R2);
#endif
}
#endif
