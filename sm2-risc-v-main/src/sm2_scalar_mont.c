/* ============================================================
 * sm2_scalar_mont.c - Montgomery-domain scalar multiplication
 *
 * Performance-test version:
 * - Public function names are unchanged.
 * - Fixed-base [k]G uses a smaller CT table: 52 windows, 31 entries/window.
 * - Scalar loops use a scalar-specific CT mixed-add to avoid the unconditional
 *   Rdbl candidate in the generic sm2_add_ja_mont().
 *
 * Important:
 * sm2_add_ja_mont_scalar_ct() is not a generic curve addition. It intentionally
 * does not cover p==q / p==-q with the generic Rdbl/Rinf candidates. Use it only
 * inside scalar multiplication paths where this risk is acceptable for testing.
 * ============================================================ */
#include "sm2_scalar.h"
#include "sm2_curve.h"
#include <string.h>

#define FIXEDBASE_W        5
#define FIXEDBASE_WINDOWS  ((256 + FIXEDBASE_W - 1) / FIXEDBASE_W)
#define FIXEDBASE_TABLE_SZ ((1 << FIXEDBASE_W) - 1)

static sm2_affine_t T_base_mont[FIXEDBASE_WINDOWS][FIXEDBASE_TABLE_SZ];
static int          T_base_mont_ready = 0;

/* ============================================================
 * Constant-time helpers
 * ============================================================ */

static int k_is_zero(const uint8_t k[32])
{
    uint8_t acc = 0;
    for (int i = 0; i < 32; i++) acc |= k[i];
    return acc == 0;
}

static inline uint32_t ct_is_zero_u32(uint32_t x)
{
    return (uint32_t)((~x & (x - 1)) >> 31);
}

static inline uint32_t ct_eq_u32(uint32_t a, uint32_t b)
{
    return ct_is_zero_u32(a ^ b);
}

static void fp_cmov(fp_t* dst, const fp_t* src, uint32_t move)
{
    uint64_t mask = 0ULL - (uint64_t)(move & 1U);
    for (int i = 0; i < 4; i++) {
        dst->v[i] = (dst->v[i] & ~mask) | (src->v[i] & mask);
    }
}

static void sm2_affine_set_infinity_local(sm2_affine_t* p)
{
    fp_set_zero(&p->x);
    fp_set_zero(&p->y);
    p->infinity = 1;
}

static void sm2_affine_cmov(sm2_affine_t* dst, const sm2_affine_t* src, uint32_t move)
{
    fp_cmov(&dst->x, &src->x, move);
    fp_cmov(&dst->y, &src->y, move);
    dst->infinity =
        (dst->infinity & (int)(~(move & 1U))) |
        (src->infinity & (int)(move & 1U));
}

void sm2_jacobian_cmov(sm2_jacobian_t* dst, const sm2_jacobian_t* src, uint32_t move)
{
    fp_cmov(&dst->X, &src->X, move);
    fp_cmov(&dst->Y, &src->Y, move);
    fp_cmov(&dst->Z, &src->Z, move);
}

static void affine_cond_neg_mont(sm2_affine_t* p, uint32_t neg)
{
    fp_t yneg;
    fp_mont_neg(&yneg, &p->y);
    fp_cmov(&p->y, &yneg, neg);
}

/* ============================================================
 * Scalar parsing helpers
 * ============================================================ */

static void scalar_be_to_u64le(uint64_t out[4], const uint8_t k[32])
{
    for (int i = 0; i < 4; i++) {
        uint64_t w = 0;
        for (int j = 0; j < 8; j++) {
            w = (w << 8) | k[i * 8 + j];
        }
        out[3 - i] = w;
    }
}

static inline uint32_t scalar_get_bits_u64le(const uint64_t a[4],
                                             unsigned bit,
                                             unsigned width)
{
    uint64_t mask = ((uint64_t)1 << width) - 1;
    unsigned limb = bit >> 6;
    unsigned shift = bit & 63;
    uint64_t v = 0;

    if (limb < 4) v = a[limb] >> shift;
    if (shift != 0 && limb + 1 < 4) v |= a[limb + 1] << (64 - shift);

    return (uint32_t)(v & mask);
}

static inline int booth_get_digit_u64le(const uint64_t a[4], unsigned w, int i)
{
    uint64_t mask = ((uint64_t)1 << w) - 1;
    uint64_t wbits;
    int n, j;

    if (i == 0) return (int)(((a[0] << 1) & mask) - (a[0] & mask));

    j = i * (int)w - 1;
    n = j / 64;
    j = j % 64;

    wbits = a[n] >> j;
    if ((64 - j) < (int)(w + 1) && n < 3) wbits |= a[n + 1] << (64 - j);

    return (int)((wbits & mask) - ((wbits >> 1) & mask));
}

/* ============================================================
 * Scalar-specific CT mixed-add
 * ============================================================ */

/*
 * r = p + q, where p is Jacobian Montgomery and q is affine Montgomery.
 *
 * This is a scalar-loop-only fast path. It keeps CT handling for p=inf and
 * q=inf, but does not build the generic Rdbl/Rinf candidates for p==q/p==-q.
 * That removes one full point doubling from every scalar-loop mixed-add.
 */
static void sm2_add_ja_mont_scalar_ct(sm2_jacobian_t* r,
                                      const sm2_jacobian_t* p,
                                      const sm2_affine_t* q_mont)
{
    uint64_t z_or;
    uint32_t pinf;
    uint32_t qinf;

    fp_t Z1Z1, Z1Z1Z1;
    fp_t U2, S2;
    fp_t H, Rv;
    fp_t HH, HHH, X1HH, X1HH2;
    fp_t RR;
    fp_t X3, Y3, Z3;
    fp_t t;

    sm2_curve_init_once_mont();

    z_or = p->Z.v[0] | p->Z.v[1] | p->Z.v[2] | p->Z.v[3];
    pinf = (uint32_t)(z_or == 0);
    qinf = (uint32_t)(q_mont->infinity != 0);

    fp_mont_sqr(&Z1Z1, &p->Z);
    fp_mont_mul(&U2, &q_mont->x, &Z1Z1);

    fp_mont_mul(&Z1Z1Z1, &Z1Z1, &p->Z);
    fp_mont_mul(&S2, &q_mont->y, &Z1Z1Z1);

    fp_mont_sub(&H,  &U2, &p->X);
    fp_mont_sub(&Rv, &S2, &p->Y);

    fp_mont_mul(&Z3, &p->Z, &H);

    fp_mont_sqr(&HH, &H);
    fp_mont_mul(&HHH, &HH, &H);
    fp_mont_mul(&X1HH, &p->X, &HH);
    fp_mont_add(&X1HH2, &X1HH, &X1HH);

    fp_mont_sqr(&RR, &Rv);

    fp_mont_sub(&t,  &RR, &HHH);
    fp_mont_sub(&X3, &t,  &X1HH2);

    fp_mont_sub(&t,  &X1HH, &X3);
    fp_mont_mul(&t,  &Rv,   &t);
    fp_mont_mul(&Y3, &p->Y, &HHH);
    fp_mont_sub(&Y3, &t,    &Y3);

    r->X = X3;
    r->Y = Y3;
    r->Z = Z3;

    /* p=inf => result is q. Inline affine_mont_to_jacobian_mont(). */
    {
        sm2_jacobian_t Rq;
        fp_t zero;

        fp_copy(&Rq.X, &q_mont->x);
        fp_copy(&Rq.Y, &q_mont->y);
        fp_copy(&Rq.Z, &FP_MONT_ONE);

        fp_set_zero(&zero);
        fp_cmov(&Rq.Z, &zero, qinf);

        sm2_jacobian_cmov(r, &Rq, pinf);
    }

    /* q=inf => result is p. */
    sm2_jacobian_cmov(r, p, qinf);
}

/* ============================================================
 * Fixed-base [k]G precomputation and multiplication
 * ============================================================ */

static void sm2_init_fixedbase_table_mont(void)
{
    if (T_base_mont_ready) return;

    sm2_affine_t base;
    sm2_get_base_affine_mont(&base);

    for (int i = 0; i < FIXEDBASE_WINDOWS; i++) {
        sm2_jacobian_t tmp[FIXEDBASE_TABLE_SZ];
        fp_t prefix_buf[FIXEDBASE_TABLE_SZ + 1];
        sm2_jacobian_t J;

        sm2_affine_mont_to_jacobian_mont(&J, &base);
        tmp[0] = J;

        for (int j = 1; j < FIXEDBASE_TABLE_SZ; j++) {
            /*
             * j==1 computes base + base, so the generic add is required.
             * j>=2 computes (j*base)+base and can use the scalar fast path.
             */
            if (j == 1) {
                sm2_add_ja_mont(&tmp[j], &tmp[j - 1], &base);
            } else {
                sm2_add_ja_mont_scalar_ct(&tmp[j], &tmp[j - 1], &base);
            }
        }

        sm2_batch_normalize_mont(T_base_mont[i], tmp, FIXEDBASE_TABLE_SZ, prefix_buf);
        for (int j = 0; j < FIXEDBASE_TABLE_SZ; j++) {
            T_base_mont[i][j].infinity = 0;
        }

        if (i + 1 < FIXEDBASE_WINDOWS) {
            sm2_jacobian_t next_base_jac;
            sm2_affine_mont_to_jacobian_mont(&next_base_jac, &base);
            for (int d = 0; d < FIXEDBASE_W; d++) {
                sm2_double_jm_a_minus3_mont(&next_base_jac, &next_base_jac);
            }
            sm2_jacobian_mont_to_affine_mont(&base, &next_base_jac);
            base.infinity = 0;
        }
    }

    T_base_mont_ready = 1;
}

static void sm2_fixedbase_select_affine_mont(sm2_affine_t* out, int window, uint32_t digit)
{
    sm2_affine_set_infinity_local(out);

    for (uint32_t v = 1; v <= FIXEDBASE_TABLE_SZ; v++) {
        uint32_t match = ct_eq_u32(digit, v);
        sm2_affine_cmov(out, &T_base_mont[window][v - 1], match);
    }
}

static void sm2_scalar_mul_fixedbase_mont(sm2_jacobian_t *R, const uint8_t k[32])
{
    uint64_t a[4];

    sm2_init_fixedbase_table_mont();
    scalar_be_to_u64le(a, k);
    sm2_jacobian_set_infinity_mont(R);

    for (int i = 0; i < FIXEDBASE_WINDOWS; i++) {
        uint32_t digit = scalar_get_bits_u64le(a, (unsigned)(i * FIXEDBASE_W), FIXEDBASE_W);
        sm2_affine_t Q;
        sm2_jacobian_t Radd;

        sm2_fixedbase_select_affine_mont(&Q, i, digit);
        sm2_add_ja_mont_scalar_ct(&Radd, R, &Q);
        *R = Radd;
    }
}

/* ============================================================
 * Generic point multiplication precomputation
 * ============================================================ */

static void sm2_precompute_table_1_to_16_affine_mont(sm2_affine_t T16[16],
                                                     const sm2_affine_t* P_mont)
{
    sm2_jacobian_t tmp[16];
    sm2_jacobian_t J;
    fp_t prefix_buf[17];

    sm2_affine_mont_to_jacobian_mont(&J, P_mont);
    tmp[0] = J;

    for (int i = 1; i < 16; i++) {
        /*
         * i==1 computes P+P, so use generic add.
         * i>=2 computes (iP)+P and can use the scalar fast path.
         */
        if (i == 1) {
            sm2_add_ja_mont(&tmp[i], &tmp[i - 1], P_mont);
        } else {
            sm2_add_ja_mont_scalar_ct(&tmp[i], &tmp[i - 1], P_mont);
        }
    }

    sm2_batch_normalize_mont(T16, tmp, 16, prefix_buf);
    for (int i = 0; i < 16; i++) {
        T16[i].infinity = 0;
    }
}

static void sm2_select_1_to_16_affine_mont(sm2_affine_t* out,
                                           const sm2_affine_t table[16],
                                           uint32_t absd)
{
    sm2_affine_set_infinity_local(out);

    for (uint32_t v = 1; v <= 16; v++) {
        uint32_t match = ct_eq_u32(absd, v);
        sm2_affine_cmov(out, &table[v - 1], match);
    }
}

/* ============================================================
 * Montgomery SchemeB core
 * ============================================================ */

static void sm2_scalar_mul_window_ct_schemeB_mont_core(
    sm2_jacobian_t* R, const sm2_affine_t* P_mont, const uint8_t k[32])
{
    if (sm2_affine_is_infinity(P_mont) || k_is_zero(k)) {
        sm2_jacobian_set_infinity_mont(R);
        return;
    }

    {
        sm2_affine_t Gm;
        uint64_t diff = 0;

        sm2_get_base_affine_mont(&Gm);
        for (int i = 0; i < 4; i++) {
            diff |= (P_mont->x.v[i] ^ Gm.x.v[i]);
            diff |= (P_mont->y.v[i] ^ Gm.y.v[i]);
        }

        if (diff == 0) {
            sm2_scalar_mul_fixedbase_mont(R, k);
            return;
        }
    }

    {
        const unsigned w = 5;
        const int n = (256 + 4) / 5;
        sm2_affine_t T16[16];
        uint64_t a[4];

        sm2_precompute_table_1_to_16_affine_mont(T16, P_mont);
        scalar_be_to_u64le(a, k);
        sm2_jacobian_set_infinity_mont(R);

        for (int i = n - 1; i >= 0; i--) {
            int di;
            int32_t d32;
            uint32_t sign;
            uint32_t absd;
            sm2_affine_t Q;
            sm2_jacobian_t Radd;

            sm2_double_jm_a_minus3_mont(R, R);
            sm2_double_jm_a_minus3_mont(R, R);
            sm2_double_jm_a_minus3_mont(R, R);
            sm2_double_jm_a_minus3_mont(R, R);
            sm2_double_jm_a_minus3_mont(R, R);

            di   = booth_get_digit_u64le(a, w, i);
            d32  = (int32_t)di;
            sign = (uint32_t)(d32 >> 31);
            absd = (uint32_t)((d32 ^ (int32_t)sign) - (int32_t)sign);

            sm2_select_1_to_16_affine_mont(&Q, T16, absd);
            affine_cond_neg_mont(&Q, sign);

            sm2_add_ja_mont_scalar_ct(&Radd, R, &Q);
            *R = Radd;
        }
    }
}

void sm2_scalar_mul_window_ct_schemeB_mont(
    sm2_jacobian_t* R, const sm2_affine_t* P, const uint8_t k[32])
{
    if (sm2_affine_is_infinity(P) || k_is_zero(k)) {
        sm2_jacobian_set_infinity_mont(R);
        return;
    }

    sm2_affine_t Pm;
    sm2_affine_to_mont(&Pm, P);
    sm2_scalar_mul_window_ct_schemeB_mont_core(R, &Pm, k);
}

void sm2_scalar_mul_window_ct_schemeB_mont_to_affine(
    sm2_affine_t* R, const sm2_affine_t* P, const uint8_t k[32])
{
    sm2_jacobian_t T;
    sm2_affine_t Pm;

    if (sm2_affine_is_infinity(P) || k_is_zero(k)) {
        R->infinity = 1;
        return;
    }

    sm2_affine_to_mont(&Pm, P);
    sm2_scalar_mul_window_ct_schemeB_mont_core(&T, &Pm, k);
    sm2_jacobian_to_affine_mont(R, &T);
}

void sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(
    sm2_jacobian_t* R, const sm2_affine_t* P_mont, const uint8_t k[32])
{
    sm2_scalar_mul_window_ct_schemeB_mont_core(R, P_mont, k);
}

void sm2_scalar_mul_window_ct_schemeB_mont_jacobian_to_jacobian_mont(
    sm2_jacobian_t* R, const sm2_jacobian_t* P_mont_jac, const uint8_t k[32])
{
    if (sm2_jacobian_is_infinity(P_mont_jac) || k_is_zero(k)) {
        sm2_jacobian_set_infinity_mont(R);
        return;
    }

    sm2_affine_t P_mont_aff;
    sm2_jacobian_mont_to_affine_mont(&P_mont_aff, P_mont_jac);
    sm2_scalar_mul_window_ct_schemeB_mont_core(R, &P_mont_aff, k);
}
