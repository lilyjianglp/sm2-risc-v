/* ============================================================
 * sm2_scalar.c 
 * ============================================================ */
#include "sm2_scalar.h"
#include <string.h>

/* ============================================================
 * scalar utils
 * ============================================================ */
static int k_is_zero(const uint8_t k[32]) {
    uint8_t acc = 0;
    for (int i = 0; i < 32; i++) acc |= k[i];
    return acc == 0;
}
/* ===== constant-time utils ===== */
static inline uint32_t ct_is_zero_u32(uint32_t x) {
    /* returns 1 if x==0 else 0 */
    return (uint32_t)((~x & (x - 1)) >> 31);
}

static inline uint32_t ct_eq_u32(uint32_t a, uint32_t b) {
    uint32_t x = a ^ b;
    return ct_is_zero_u32(x);
}

static inline uint32_t ct_is_nonzero_u32(uint32_t x) {
    return 1U ^ ct_is_zero_u32(x);
}

/* conditional move for fp_t: if move==1, dst=src else keep dst */
static void fp_cmov(fp_t* dst, const fp_t* src, uint32_t move)
{
    uint64_t mask = 0ULL - (uint64_t)(move & 1U);
    for (int i = 0; i < 4; i++) {
        uint64_t d = dst->v[i];
        uint64_t s = src->v[i];
        dst->v[i] = (d & ~mask) | (s & mask);
    }
}

static void sm2_affine_cmov(sm2_affine_t* dst, const sm2_affine_t* src, uint32_t move)
{
    fp_cmov(&dst->x, &src->x, move);
    fp_cmov(&dst->y, &src->y, move);
    /* infinity 也要 cmov，但这里保证表项都不是 infinity */
    dst->infinity = (dst->infinity & (int)(~(move & 1U))) | (src->infinity & (int)(move & 1U));
}

void sm2_jacobian_cmov(sm2_jacobian_t* dst, const sm2_jacobian_t* src, uint32_t move)
{
    fp_cmov(&dst->X, &src->X, move);
    fp_cmov(&dst->Y, &src->Y, move);
    fp_cmov(&dst->Z, &src->Z, move);
}

/* conditional negation on affine y: if neg==1 => y = -y */
static void affine_cond_neg(sm2_affine_t* p, uint32_t neg)
{
    fp_t yneg;
    fp_neg(&yneg, &p->y);
    fp_cmov(&p->y, &yneg, neg);
}
static void scalar_be_to_u64le(uint64_t out[4], const uint8_t k[32])
{
    /* k[0] is MSB, k[31] is LSB */
    for (int i = 0; i < 4; i++) {
        uint64_t w = 0;
        for (int j = 0; j < 8; j++) {
            w = (w << 8) | k[i*8 + j];
        }
        out[3 - i] = w; /* little-endian limbs: out[0] is low */
    }
}

static inline int booth_get_digit_u64le(const uint64_t a[4], unsigned w, int i)
{
    /* same math as gmssl get_booth, i is public */
    uint64_t mask = ((uint64_t)1 << w) - 1;//低位连续有w个1
    uint64_t wbits;//w+1位
    int n, j;

    if (i == 0) {
        return (int)(((a[0] << 1) & mask) - (a[0] & mask));
    }

    j = i * (int)w - 1;
    n = j /64;
    j = j %64;

    wbits = a[n] >> j;
    if ((64 - j) < (int)(w + 1) && n < 3) {
        wbits |= a[n + 1] << (64 - j);
    }

    return (int)((wbits & mask) - ((wbits >> 1) & mask));
}


/* build table: T[i] = (i+1)*P, i=0..15  (Affine) */
static void sm2_precompute_table_1_to_16_affine(sm2_affine_t T16[16],
                                                const sm2_affine_t* P)
{
    sm2_jacobian_t tmp[16];
    sm2_jacobian_t J;

    sm2_affine_to_jacobian(&J, P);

    /* tmp[0] = 1P */
    tmp[0] = J;

    /* tmp[i] = tmp[i-1] + P  (use mixed add JA) */
    for (int i = 1; i < 16; i++) {
        sm2_add_ja(&tmp[i], &tmp[i-1], P);
    }

    /* batch normalize to affine */
    fp_t prefix_buf[16 + 1];
    sm2_batch_normalize(T16, tmp, 16, prefix_buf);


    for (int i = 0; i < 16; i++) T16[i].infinity = 0;
}

/* ============================================================
 * constant-time helpers (ONLY used by ladder)
 * ============================================================ */
static inline uint32_t scalar_get_bit_be(const uint8_t k[32], int bit)
{
    /* bit: 255 = MSB, 0 = LSB */
    int byte = bit >> 3;
    int off = bit & 7;
    return (k[31 - byte] >> off) & 1U;
}

static void fp_cswap(fp_t* a, fp_t* b, uint32_t swap)
{
    uint64_t mask = (uint64_t)0 - (uint64_t)(swap & 1U);
    for (int i = 0; i < 4; i++) {
        uint64_t t = (a->v[i] ^ b->v[i]) & mask;
        a->v[i] ^= t;
        b->v[i] ^= t;
    }
}

void sm2_jacobian_cswap(sm2_jacobian_t* a,
    sm2_jacobian_t* b,
    uint32_t swap)
{
    fp_cswap(&a->X, &b->X, swap);
    fp_cswap(&a->Y, &b->Y, swap);
    fp_cswap(&a->Z, &b->Z, swap);
}

/* ============================================================
 * 1) NAF scalar multiplication (non-constant-time)����
 * ============================================================ */
static int k_is_odd(const uint8_t k[32]) {
    return (k[31] & 1) != 0;
}

static int k_mod4(const uint8_t k[32]) {
    return (int)(k[31] & 3);
}

static void k_add1(uint8_t k[32]) {
    for (int i = 31; i >= 0; i--) {
        if (++k[i] != 0) break;
    }
}

static void k_sub1(uint8_t k[32]) {
    for (int i = 31; i >= 0; i--) {
        if (k[i]-- != 0) break;
    }
}

static void k_rshift1(uint8_t k[32]) {
    uint8_t carry = 0;
    for (int i = 0; i < 32; i++) {
        uint8_t next = (k[i] & 1) << 7;
        k[i] = (k[i] >> 1) | carry;
        carry = next;
    }
}
/*w=2*/
static size_t naf_encode_i8(int8_t naf[257], const uint8_t k_be[32])
{
    uint8_t k[32];
    memcpy(k, k_be, 32);

    size_t i = 0;
    while (!k_is_zero(k)) {
        if (!k_is_odd(k)) {
            naf[i++] = 0;
        }
        else {
            int8_t di = (int8_t)(2 - k_mod4(k)); /* +1 or -1 */
            naf[i++] = di;
            if (di == 1) k_sub1(k);
            else         k_add1(k);
        }
        k_rshift1(k);
    }
    return i;
}

void sm2_scalar_mul_naf(sm2_jacobian_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32])
{
    int8_t naf[257];
    size_t len = naf_encode_i8(naf, k);

    sm2_jacobian_set_infinity(R);
    if (len == 0 || sm2_affine_is_infinity(P)) return;

    sm2_affine_t Pneg = *P;
    fp_neg(&Pneg.y, &P->y);

    for (size_t i = len; i-- > 0;) {
        sm2_double_jm_a_minus3(R, R);
        if (naf[i] == 1)
            sm2_add_ja(R, R, P);
        else if (naf[i] == -1)
            sm2_add_ja(R, R, &Pneg);
    }
}

/* ============================================================
 * 2) Constant-time Montgomery ladder
 * ============================================================ */
void sm2_scalar_mul_ladder_ct(sm2_jacobian_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32])
{
    sm2_jacobian_t R0, R1, T;

    if (sm2_affine_is_infinity(P) || k_is_zero(k)) {
        sm2_jacobian_set_infinity(R);
        return;
    }

    sm2_jacobian_set_infinity(&R0);
    sm2_affine_to_jacobian(&R1, P);

    uint32_t prev = 0;
    for (int i = 255; i >= 0; i--) {
        uint32_t b = scalar_get_bit_be(k, i);
        uint32_t swap = b ^ prev;
        sm2_jacobian_cswap(&R0, &R1, swap);
        prev = b;

        sm2_add_jj(&T, &R0, &R1);
        sm2_double_jm_a_minus3(&R0, &R0);
        R1 = T;
    }

    sm2_jacobian_cswap(&R0, &R1, prev);
    *R = R0;
}
/* ============================================================
 *   SchemeB: fixed w=5, fully-unrolled 5 doubles, CT select
 *    Works for BOTH odd/even k now.
 * ============================================================ */
void sm2_scalar_mul_window_ct_schemeB(sm2_jacobian_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32])
{
    const unsigned w = 5;
    const int n = (256 + 4) / 5; /* 52 */

    if (sm2_affine_is_infinity(P) || k_is_zero(k)) {
        sm2_jacobian_set_infinity(R);
        return;
    }

    /* 1) 预计算 1..16 表（Affine） */
    sm2_affine_t T16[16];
    sm2_precompute_table_1_to_16_affine(T16, P);

    /* 2) scalar -> u64 little-endian limbs */
    uint64_t a[4];
    scalar_be_to_u64le(a, k);

    sm2_jacobian_set_infinity(R);

    for (int i = n - 1; i >= 0; i--) {

        /* ===== fully unrolled: 5 doubles ===== */
        sm2_double_jm_a_minus3(R, R);
        sm2_double_jm_a_minus3(R, R);
        sm2_double_jm_a_minus3(R, R);
        sm2_double_jm_a_minus3(R, R);
        sm2_double_jm_a_minus3(R, R);

        /* ===== Booth digit (may be even) ===== */
        int di =booth_get_digit_u64le(a, w, i);
;

        int32_t d32 = (int32_t)di;
        uint32_t sign = (uint32_t)(d32 >> 31);
        uint32_t absd = (uint32_t)((d32 ^ (int32_t)sign) - (int32_t)sign); /* abs */
        uint32_t nz   = ct_is_nonzero_u32(absd);

        /* ===== constant-time select Q = absd*P from T16 =====
           absd in {0..15}. absd==0 => later cmov keeps R.
         */
        sm2_affine_t Q = T16[0];
        Q.infinity = 0;

        /* select absd == 1..15 => T16[absd-1] */
        for (uint32_t v = 1; v <= 16; v++) {
            uint32_t match = ct_eq_u32(absd, v);
            sm2_affine_cmov(&Q, &T16[v - 1], match);
        }

        /* sign fix */
        affine_cond_neg(&Q, sign);

        /* Radd = R + Q */
        sm2_jacobian_t Radd;
        sm2_add_ja(&Radd, R, &Q);

        /* if absd==0 keep R else take Radd */
        sm2_jacobian_cmov(R, &Radd, nz);
    }
}
