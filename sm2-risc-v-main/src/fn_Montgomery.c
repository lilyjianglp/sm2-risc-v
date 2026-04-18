#include "fn.h"
#include <string.h>
#include <stdint.h>
/* ============================================================
 *  Constants: SM2 curve order n
 *  n = FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123
 * ============================================================ */
static inline void fn_cond_sub_n(fn_t *r);

static const fn_t FN_N = {
    .v = {
        0x53BBF40939D54123ULL,
        0x7203DF6B21C6052BULL,
        0xFFFFFFFFFFFFFFFFULL,
        0xFFFFFFFEFFFFFFFFULL
    }
};


const fn_t* fn_modulus_n(void) {
    return &FN_N;
}

/* ============================================================
 *  Basic helpers
 * ============================================================ */

void fn_zero(fn_t *r) {
    r->v[0] = r->v[1] = r->v[2] = r->v[3] = 0;
}

int fn_is_zero(const fn_t *a) {
    return (a->v[0] | a->v[1] | a->v[2] | a->v[3]) == 0;
}

int fn_cmp(const fn_t *a, const fn_t *b) {
    for (int i = 3; i >= 0; i--) {
        if (a->v[i] < b->v[i]) return -1;
        if (a->v[i] > b->v[i]) return 1;
    }
    return 0;
}


static inline uint64_t addc_u64(uint64_t a, uint64_t b, uint64_t *carry) {
    uint64_t c = *carry;
    uint64_t r = a + b;
    uint64_t c1 = (r < a);         // carry from a+b
    r = r + c;
    uint64_t c2 = (r < c);         // carry from +c
    *carry = c1 | c2;
    return r;
}

static inline uint64_t subb_u64(uint64_t a, uint64_t b, uint64_t *borrow) {
    uint64_t c = *borrow;
    uint64_t r = a - b;
    uint64_t b1 = (a < b);         // borrow from a-b
    uint64_t r2 = r - c;
    uint64_t b2 = (r < c);         // borrow from (a-b)-c
    *borrow = b1 | b2;
    return r2;
}



void fn_from_be(fn_t *r, const uint8_t in[32]) {
    for (int i = 0; i < 4; i++) {
        uint64_t w = 0;
        for (int j = 0; j < 8; j++) {
            w = (w << 8) | (uint64_t)in[i * 8 + j];
        }
        r->v[3 - i] = w;
    }
    /* 关键：一次规约，保证 < n */
    fn_cond_sub_n(r);
}


void fn_to_be(uint8_t out[32], const fn_t *a) {
    for (int i = 0; i < 4; i++) {
        uint64_t w = a->v[3 - i];
        out[i * 8 + 0] = (uint8_t)(w >> 56);
        out[i * 8 + 1] = (uint8_t)(w >> 48);
        out[i * 8 + 2] = (uint8_t)(w >> 40);
        out[i * 8 + 3] = (uint8_t)(w >> 32);
        out[i * 8 + 4] = (uint8_t)(w >> 24);
        out[i * 8 + 5] = (uint8_t)(w >> 16);
        out[i * 8 + 6] = (uint8_t)(w >> 8);
        out[i * 8 + 7] = (uint8_t)(w);
    }
}






/* ============================================================
 *  Arithmetic mod n
 * ============================================================ */

static inline uint64_t ct_mask_u64(uint64_t bit) {
    return (uint64_t)0 - (bit & 1ULL);
}

static const fn_t FN_C = { /* C = 2^256 - n */
    .v = {
        0xAC440BF6C62ABEDDULL,
        0x8DFC2094DE39FAD4ULL,
        0x0000000000000000ULL,
        0x0000000100000000ULL
    }
};

void fn_add(fn_t *r, const fn_t *a, const fn_t *b) {
    fn_t t;

    uint64_t carry = 0;
    t.v[0] = addc_u64(a->v[0], b->v[0], &carry);
    t.v[1] = addc_u64(a->v[1], b->v[1], &carry);
    t.v[2] = addc_u64(a->v[2], b->v[2], &carry);
    t.v[3] = addc_u64(a->v[3], b->v[3], &carry);

    /* branchless: if carry then t += C */
    uint64_t mcarry = ct_mask_u64(carry);
    uint64_t c2 = 0;
    t.v[0] = addc_u64(t.v[0], FN_C.v[0] & mcarry, &c2);
    t.v[1] = addc_u64(t.v[1], FN_C.v[1] & mcarry, &c2);
    t.v[2] = addc_u64(t.v[2], FN_C.v[2] & mcarry, &c2);
    t.v[3] = addc_u64(t.v[3], FN_C.v[3] & mcarry, &c2);

    /* 再做一次条件减 n（固定次数） */
    fn_cond_sub_n(&t);
    *r = t;
}



void fn_sub(fn_t *r, const fn_t *a, const fn_t *b) {
    fn_t t;
    uint64_t borrow = 0;

    t.v[0] = subb_u64(a->v[0], b->v[0], &borrow);
    t.v[1] = subb_u64(a->v[1], b->v[1], &borrow);
    t.v[2] = subb_u64(a->v[2], b->v[2], &borrow);
    t.v[3] = subb_u64(a->v[3], b->v[3], &borrow);

    /* if borrow==1 => add n back (branchless) */
    uint64_t m = ct_mask_u64(borrow);
    uint64_t carry = 0;
    t.v[0] = addc_u64(t.v[0], FN_N.v[0] & m, &carry);
    t.v[1] = addc_u64(t.v[1], FN_N.v[1] & m, &carry);
    t.v[2] = addc_u64(t.v[2], FN_N.v[2] & m, &carry);
    t.v[3] = addc_u64(t.v[3], FN_N.v[3] & m, &carry);

    *r = t; /* 已经在 [0,n) */
}


/* ============================================================
 * Montgomery arithmetic mod n (SM2 order)
 * radix b = 2^64, k = 4 limbs, R = b^k = 2^256
 * ============================================================ */

/* ============================================================
 * Constant-time Montgomery arithmetic mod n (SM2 order)
 * radix b = 2^64, k = 4 limbs, R = b^k = 2^256
 * ============================================================ */

static const uint64_t FN_N0PRIME = 0x327f9e8872350975ULL; /* -n^{-1} mod 2^64 */

/* mont(1) = R mod n = 2^256 - n */
static const fn_t FN_MONT_ONE = {
    .v = {
        0xAC440BF6C62ABEDDULL,
        0x8DFC2094DE39FAD4ULL,
        0x0000000000000000ULL,
        0x0000000100000000ULL
    }
};

/* R^2 mod n = 2^512 mod n */
static const fn_t FN_MONT_R2 = {
    .v = {
        0x901192AF7C114F20ULL,
        0x3464504ADE6FA2FAULL,
        0x620FC84C3AFFE0D4ULL,
        0x1EB5E412A22B3D3BULL
    }
};

/* r = r - n if r >= n  (branchless) */
static inline void fn_cond_sub_n(fn_t *r) {
    uint64_t br = 0;
    fn_t t;

    t.v[0] = subb_u64(r->v[0], FN_N.v[0], &br);
    t.v[1] = subb_u64(r->v[1], FN_N.v[1], &br);
    t.v[2] = subb_u64(r->v[2], FN_N.v[2], &br);
    t.v[3] = subb_u64(r->v[3], FN_N.v[3], &br);

    /* br==0 => r>=n => take t */
    uint64_t mask = ct_mask_u64(br ^ 1ULL);
    r->v[0] = (t.v[0] & mask) | (r->v[0] & ~mask);
    r->v[1] = (t.v[1] & mask) | (r->v[1] & ~mask);
    r->v[2] = (t.v[2] & mask) | (r->v[2] & ~mask);
    r->v[3] = (t.v[3] & mask) | (r->v[3] & ~mask);
}

/* if (t8==1) r += (R mod n)  (branchless) */
static inline void fn_cond_add_R_if_t8(fn_t *r, uint64_t t8) {
    uint64_t mask = ct_mask_u64(t8 & 1ULL);
    uint64_t carry = 0;

    r->v[0] = addc_u64(r->v[0], FN_MONT_ONE.v[0] & mask, &carry);
    r->v[1] = addc_u64(r->v[1], FN_MONT_ONE.v[1] & mask, &carry);
    r->v[2] = addc_u64(r->v[2], FN_MONT_ONE.v[2] & mask, &carry);
    r->v[3] = addc_u64(r->v[3], FN_MONT_ONE.v[3] & mask, &carry);
    (void)carry;
}




/* 4x4 -> 8 limbs multiply (512-bit), constant-time (fixed carry propagation) */


/* 128-bit (hi:lo) += x (x is 64-bit) */
static inline void add_u64_to_u128(uint64_t *lo, uint64_t *hi, uint64_t x)
{
    uint64_t c = 0;
    *lo = addc_u64(*lo, x, &c);
    *hi = *hi + c;
}

/*
 * 64x64 -> 128 multiply without __int128, using 32-bit decomposition.
 *
 * Let a = a1*2^32 + a0, b = b1*2^32 + b0
 * Then:
 *   a*b = (a0*b0)
 *       + (a0*b1 + a1*b0) * 2^32
 *       + (a1*b1) * 2^64
 *
 * We compute lo and hi in base 2^32 with carries.
 */
 
 

static inline void mul_u64_to_u128_ref(uint64_t a, uint64_t b,
                                           uint64_t *lo, uint64_t *hi)
{
    const uint64_t MASK32 = 0xFFFFFFFFULL;

    uint64_t a0 = a & MASK32;
    uint64_t a1 = a >> 32;
    uint64_t b0 = b & MASK32;
    uint64_t b1 = b >> 32;

    uint64_t p00 = a0 * b0;  /* < 2^64 */
    uint64_t p01 = a0 * b1;  /* < 2^64 */
    uint64_t p10 = a1 * b0;  /* < 2^64 */
    uint64_t p11 = a1 * b1;  /* < 2^64 */

    uint64_t p00_lo = p00 & MASK32;
    uint64_t p00_hi = p00 >> 32;

    uint64_t mid = p00_hi + (p01 & MASK32) + (p10 & MASK32);
    uint64_t mid_lo = mid & MASK32;
    uint64_t mid_hi = mid >> 32;

    *lo = p00_lo | (mid_lo << 32);
    *hi = p11 + (p01 >> 32) + (p10 >> 32) + mid_hi;
}


static inline void mul_u64_to_u128(uint64_t a, uint64_t b,
                                   uint64_t *lo, uint64_t *hi)
{
    /* C 参考实现：无 __int128 */
    mul_u64_to_u128_ref(a, b, lo, hi);
}



/*
 * Equivalent to your original fn_mul_4x4_u512_ct but without __int128.
 * out[0..7] = a(4 limbs) * b(4 limbs), limbs are uint64 little-endian.
 */
#ifndef USE_ASM_FN_MONT
  __attribute__((weak))
void fn_mul_4x4_u512_ct(uint64_t out[8], const fn_t *a, const fn_t *b)
{
    for (int i = 0; i < 8; i++) out[i] = 0;

    for (int i = 0; i < 4; i++) {
        uint64_t carry = 0;

        for (int j = 0; j < 4; j++) {
            uint64_t lo, hi;

            /* (hi:lo) = a[i]*b[j] */
            mul_u64_to_u128(a->v[i], b->v[j], &lo, &hi);

            /* (hi:lo) += out[i+j] */
            add_u64_to_u128(&lo, &hi, out[i + j]);

            /* (hi:lo) += carry */
            add_u64_to_u128(&lo, &hi, carry);

            out[i + j] = lo;
            carry      = hi;
        }

        /* fixed propagation to out[i+4..7] */
        for (int k = i + 4; k < 8; k++) {
            uint64_t c = 0;
            out[k] = addc_u64(out[k], carry, &c);
            carry  = c; /* 0/1 */
        }
    }
}


/* REDC: t[0..8] -> r (constant-time carry propagation) */
__attribute__((weak))
void fn_mont_reduce_ct(fn_t *r, uint64_t t[9])
{
    for (int i = 0; i < 4; i++) {
        /* m = (t[i] * n0prime) mod 2^64 ：只要 low 64 */
        uint64_t m_lo, m_hi;
        mul_u64_to_u128(t[i], FN_N0PRIME, &m_lo, &m_hi);
        (void)m_hi;
        uint64_t m = m_lo;

        uint64_t carry = 0;

        for (int j = 0; j < 4; j++) {
            uint64_t lo, hi;

            /* (hi:lo) = m * n[j] */
            mul_u64_to_u128(m, FN_N.v[j], &lo, &hi);

            /* (hi:lo) += t[i+j] */
            add_u64_to_u128(&lo, &hi, t[i + j]);

            /* (hi:lo) += carry */
            add_u64_to_u128(&lo, &hi, carry);

            t[i + j] = lo;
            carry    = hi;
        }

        /* fixed propagation to t[i+4..8] */
        for (int k = i + 4; k < 9; k++) {
            uint64_t c = 0;
            t[k] = addc_u64(t[k], carry, &c);
            carry = c; /* 0/1 */
        }
    }

    r->v[0] = t[4];
    r->v[1] = t[5];
    r->v[2] = t[6];
    r->v[3] = t[7];

    fn_cond_add_R_if_t8(r, t[8]);

    fn_cond_sub_n(r);
    fn_cond_sub_n(r);
}


/* MontMul: r = a*b*R^{-1} mod n  (inputs/outputs in Montgomery domain) */
__attribute__((weak))
void fn_mont_mul_ct(fn_t *r, const fn_t *a, const fn_t *b) {
    uint64_t prod[8];
    uint64_t t[9];

    fn_mul_4x4_u512_ct(prod, a, b);
    for (int i = 0; i < 8; i++) t[i] = prod[i];
    t[8] = 0; /* IMPORTANT: only REDC may set t[8] */

    fn_mont_reduce_ct(r, t);
}

#endif

/* Convert into Montgomery domain: a -> aR mod n = MontMul(a, R^2 mod n) */
void fn_to_mont(fn_t *r, const fn_t *a) {
    fn_mont_mul_ct(r, a, &FN_MONT_R2);
}

static const fn_t FN_ONE = { .v = {1,0,0,0} };

/* Convert out: aR -> a = MontMul(aR, 1) */
void fn_from_mont(fn_t *r, const fn_t *a) {
    fn_mont_mul_ct(r, a, &FN_ONE);
}

/* Public: r = a*b mod n, constant-time-ish under canonical input contract */
void fn_mul(fn_t *r, const fn_t *a, const fn_t *b) {
    fn_t aa = *a, bb = *b;

    /* If you want strict contract: remove these two lines and require caller canonical.
       Keeping them is still branchless and fixed time (just extra work). */
    fn_cond_sub_n(&aa);
    fn_cond_sub_n(&bb);

    fn_t ma, mb, mr;
    fn_to_mont(&ma, &aa);
    fn_to_mont(&mb, &bb);
    fn_mont_mul_ct(&mr, &ma, &mb);
    fn_from_mont(r, &mr);
}




void fn_mont_one(fn_t *r) {
    *r = FN_MONT_ONE;
}

void fn_to_mont_pub(fn_t *r, const fn_t *a) {
    fn_t t = *a;
    fn_cond_sub_n(&t);   /* 规范化，固定时间 */
    fn_to_mont(r, &t);
}

void fn_from_mont_pub(fn_t *r, const fn_t *a) {
    fn_t t = *a;
    fn_cond_sub_n(&t);   /* 保守处理 */
    fn_from_mont(r, &t);
}

void fn_mont_mul_pub(fn_t *r, const fn_t *a, const fn_t *b) {
    fn_t aa = *a, bb = *b;
    fn_cond_sub_n(&aa);
    fn_cond_sub_n(&bb);
    fn_mont_mul_ct(r, &aa, &bb);
}

void fn_mont_sqr_pub(fn_t *r, const fn_t *a) {
    fn_mont_mul_pub(r, a, a);
}








/* ============================================================
 *  KEX helpers
 * ============================================================ */

void fn_x_dash_w127(fn_t *x_dash, const uint8_t x_be[32]) {
    fn_t x;
    fn_from_be(&x, x_be);

    /* keep low 127 bits: v[0] full, low 63 bits of v[1] */
    x.v[1] &= 0x7FFFFFFFFFFFFFFFULL;
    x.v[2] = 0;
    x.v[3] = 0;

    /* add 2^127 -> set bit 63 of limb1 */
    x.v[1] |= 0x8000000000000000ULL;

    *x_dash = x;
}

/* 0/1 -> 0x00..00 or 0xFF..FF */
static inline uint8_t ct_mask_u8(uint8_t bit) {
    return (uint8_t)(0u - (uint32_t)(bit & 1u));
}

/* constant-time: returns 1 if a==0 else 0 */
static inline uint8_t fn_is_zero_ct_u8(const fn_t *a) {
    /* fn_is_zero 本身是无分支的，这里只是显式转成 0/1 */
    return (uint8_t)(fn_is_zero(a) ? 1 : 0);
}

/* constant-time select: dst = (mask? a : b), mask=0xFF choose a, mask=0 choose b */
static inline void ct_select_u8_buf(uint8_t *dst,
                                    const uint8_t *a,
                                    const uint8_t *b,
                                    size_t n,
                                    uint8_t mask)
{
    for (size_t i = 0; i < n; i++) {
        uint8_t ai = a[i];
        uint8_t bi = b[i];
        dst[i] = (uint8_t)((ai & mask) | (bi & (uint8_t)~mask));
    }
}

int fn_compute_t(uint8_t t_be[32],
                    const uint8_t d_be[32],
                    const uint8_t r_be[32],
                    const fn_t *x)
{
    fn_t d, rr, xr, t;

    fn_from_be(&d, d_be);
    fn_from_be(&rr, r_be);

    /* ok = 1 initially, then AND with conditions */
    uint8_t ok = 1;

    /* ok &= (d != 0) */
    ok &= (uint8_t)(1u ^ fn_is_zero_ct_u8(&d));
    /* ok &= (r != 0) */
    ok &= (uint8_t)(1u ^ fn_is_zero_ct_u8(&rr));

    /* Always do the arithmetic (no early return) */
    fn_mul(&xr, x, &rr);
    fn_add(&t, &d, &xr);

    /* ok &= (t != 0) */
    ok &= (uint8_t)(1u ^ fn_is_zero_ct_u8(&t));

    /* Serialize t unconditionally */
    uint8_t tmp[32];
    fn_to_be(tmp, &t);

    /* If ok==1 => output tmp, else output all-zero */
    uint8_t zero[32] = {0};
    uint8_t mask = ct_mask_u8(ok);     /* 0xFF if ok else 0x00 */
    ct_select_u8_buf(t_be, tmp, zero, 32, mask);

    return (int)ok;
}


