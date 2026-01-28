#include "fn.h"
#include <string.h>

/* ============================================================
 *  Constants: SM2 curve order n
 *  n = FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123
 * ============================================================ */

static const fn_t FN_N = {
    .v = {
        0x53BBF40939D54123ULL,
        0x7203DF6B21C6052BULL,
        0xFFFFFFFFFFFFFFFFULL,
        0xFFFFFFFEFFFFFFFFULL
    }
};

/* Barrett mu = floor(b^(2k) / n), b=2^64, k=4 => 2k=8 limbs, mu is (k+1)=5 limbs.
 * This mu matches the one we used in kex: it's correct for SM2 n.
 */
typedef struct { uint64_t v[5]; } u320_t;
typedef struct { uint64_t v[8]; } u512_t;

static const u320_t FN_N_MU = {
    .v = {
        0x12AC6361F15149A0ULL,
        0x8DFC2096FA323C01ULL,
        0x0000000100000001ULL,
        0x0000000100000001ULL,
        0x0000000000000001ULL
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


static void fn_add_raw(fn_t *r, const fn_t *a, const fn_t *b) {
    uint64_t carry = 0;
    r->v[0] = addc_u64(a->v[0], b->v[0], &carry);
    r->v[1] = addc_u64(a->v[1], b->v[1], &carry);
    r->v[2] = addc_u64(a->v[2], b->v[2], &carry);
    r->v[3] = addc_u64(a->v[3], b->v[3], &carry);
    /* ignore overflow */
}

static void fn_sub_raw(fn_t *r, const fn_t *a, const fn_t *b, uint64_t *borrow_out) {
    uint64_t borrow = 0;
    r->v[0] = subb_u64(a->v[0], b->v[0], &borrow);
    r->v[1] = subb_u64(a->v[1], b->v[1], &borrow);
    r->v[2] = subb_u64(a->v[2], b->v[2], &borrow);
    r->v[3] = subb_u64(a->v[3], b->v[3], &borrow);
    if (borrow_out) *borrow_out = borrow;
}

void fn_from_be(fn_t *r, const uint8_t in[32]) {
    /* big-endian bytes -> little-endian 64-bit limbs */
    for (int i = 0; i < 4; i++) {
        uint64_t w = 0;
        for (int j = 0; j < 8; j++) {
            w = (w << 8) | (uint64_t)in[i * 8 + j];
        }
        r->v[3 - i] = w;
    }
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

void fn_normalize(fn_t *a) {
    /* ensure a in [0,n) by repeated subtraction (fast here: values are usually already small-ish) */
    while (fn_cmp(a, &FN_N) >= 0) {
        fn_t t = *a;
        uint64_t borrow = 0;
        a->v[0] = subb_u64(t.v[0], FN_N.v[0], &borrow);
        a->v[1] = subb_u64(t.v[1], FN_N.v[1], &borrow);
        a->v[2] = subb_u64(t.v[2], FN_N.v[2], &borrow);
        a->v[3] = subb_u64(t.v[3], FN_N.v[3], &borrow);
    }
}

/* ============================================================
 *  256x256 -> 512 multiply
 * ============================================================ */

static void mul_4x4(u512_t *out, const fn_t *a, const fn_t *b) {
    uint64_t r[8] = {0};

    for (int i = 0; i < 4; i++) {
        unsigned __int128 carry = 0;
        for (int j = 0; j < 4; j++) {
            unsigned __int128 t = (unsigned __int128)a->v[i] * b->v[j] + r[i + j] + carry;
            r[i + j] = (uint64_t)t;
            carry = t >> 64;
        }
        r[i + 4] += (uint64_t)carry;
    }

    for (int i = 0; i < 8; i++) out->v[i] = r[i];
}

/* ============================================================
 *  Barrett reduction: x (512-bit) -> x mod n
 *  Using b=2^64, k=4:
 *    q1 = floor(x / b^(k-1))  -> shift right 192 bits -> 5 limbs
 *    q2 = q1 * mu             -> 10 limbs
 *    q3 = floor(q2 / b^(k+1)) -> shift right 320 bits -> 5 limbs
 *    r1 = x mod b^(k+1)       -> low 5 limbs
 *    r2 = (q3 * n) mod b^(k+1)-> low 5 limbs
 *    r  = r1 - r2             -> mod b^(k+1)
 *    while r >= n: r -= n     -> final
 * ============================================================ */

static void mul_5x5_10(uint64_t out10[10], const uint64_t a[5], const uint64_t b[5]) {
    for (int i = 0; i < 10; i++) out10[i] = 0;

    for (int i = 0; i < 5; i++) {
        unsigned __int128 carry = 0;
        for (int j = 0; j < 5; j++) {
            unsigned __int128 t = (unsigned __int128)a[i] * b[j] + out10[i + j] + carry;
            out10[i + j] = (uint64_t)t;
            carry = t >> 64;
        }
        out10[i + 5] += (uint64_t)carry;
    }
}



static void barrett_mod_n(fn_t *r, const u512_t *x) {
    /* q1 = floor(x / b^(k-1)) with k=4 => shift right 3 limbs => take limbs[3..7] */
    uint64_t q1[5];
    for (int i = 0; i < 5; i++) q1[i] = x->v[i + 3];

    /* q2 = q1 * mu (10 limbs) */
    uint64_t q2[10];
    mul_5x5_10(q2, q1, FN_N_MU.v);

    /* q3 = floor(q2 / b^(k+1)) => shift right 5 limbs => take q2[5..9] */
    uint64_t q3[5];
    for (int i = 0; i < 5; i++) q3[i] = q2[i + 5];

    /* r1 = x mod b^(k+1) => low 5 limbs */
    uint64_t r1[5];
    for (int i = 0; i < 5; i++) r1[i] = x->v[i];

    /* r2 = (q3 * n) mod b^(k+1) => compute (5 limbs)*(4 limbs), keep low 5 limbs */
    uint64_t tmp[9] = {0};
    for (int i = 0; i < 5; i++) {
        unsigned __int128 carry = 0;
        for (int j = 0; j < 4; j++) {
            unsigned __int128 t = (unsigned __int128)q3[i] * FN_N.v[j] + tmp[i + j] + carry;
            tmp[i + j] = (uint64_t)t;
            carry = t >> 64;
        }
        tmp[i + 4] += (uint64_t)carry;
    }

    uint64_t r2[5];
    for (int i = 0; i < 5; i++) r2[i] = tmp[i];

    /* rr = r1 - r2 (mod 2^(320))  —— 用 subb_u64 正确处理借位链 */
    uint64_t rr[5];
    uint64_t borrow = 0;
    for (int i = 0; i < 5; i++) {
        rr[i] = subb_u64(r1[i], r2[i], &borrow);
    }
    (void)borrow; /* wrap-around is fine (mod 2^320) */

    /* Reduce rr (5 limbs) by n (4 limbs) until rr < n */
    for (;;) {
        int ge;

        if (rr[4] != 0) {
            ge = 1;
        } else {
            /* rr[3..0] >= n ? */
            if (rr[3] > FN_N.v[3]) ge = 1;
            else if (rr[3] < FN_N.v[3]) ge = 0;
            else if (rr[2] > FN_N.v[2]) ge = 1;
            else if (rr[2] < FN_N.v[2]) ge = 0;
            else if (rr[1] > FN_N.v[1]) ge = 1;
            else if (rr[1] < FN_N.v[1]) ge = 0;
            else ge = (rr[0] >= FN_N.v[0]);
        }

        if (!ge) break;

        uint64_t br = 0;
        rr[0] = subb_u64(rr[0], FN_N.v[0], &br);
        rr[1] = subb_u64(rr[1], FN_N.v[1], &br);
        rr[2] = subb_u64(rr[2], FN_N.v[2], &br);
        rr[3] = subb_u64(rr[3], FN_N.v[3], &br);
        rr[4] = subb_u64(rr[4], 0,        &br);
    }

    r->v[0] = rr[0];
    r->v[1] = rr[1];
    r->v[2] = rr[2];
    r->v[3] = rr[3];
}





/* ============================================================
 *  Arithmetic mod n
 * ============================================================ */

void fn_add(fn_t *r, const fn_t *a, const fn_t *b) {
    fn_t t;
    uint64_t carry = 0;

    t.v[0] = addc_u64(a->v[0], b->v[0], &carry);
    t.v[1] = addc_u64(a->v[1], b->v[1], &carry);
    t.v[2] = addc_u64(a->v[2], b->v[2], &carry);
    t.v[3] = addc_u64(a->v[3], b->v[3], &carry);

    /* 如果发生 256-bit 溢出：加上 C = 2^256 - n （因为 2^256 ≡ C (mod n)） */
    if (carry) {
        static const fn_t C = {
            .v = {
                0xAC440BF6C62ABEDDULL,
                0x8DFC2094DE39FAD4ULL,
                0x0000000000000000ULL,
                0x0000000100000000ULL
            }
        };
        uint64_t c2 = 0;
        t.v[0] = addc_u64(t.v[0], C.v[0], &c2);
        t.v[1] = addc_u64(t.v[1], C.v[1], &c2);
        t.v[2] = addc_u64(t.v[2], C.v[2], &c2);
        t.v[3] = addc_u64(t.v[3], C.v[3], &c2);
        /* c2 溢出可以忽略：因为 (t + C) < 2n，最终再规约即可 */
    }

    /* 最终规约到 [0,n) */
    if (fn_cmp(&t, &FN_N) >= 0) {
        uint64_t borrow = 0;
        t.v[0] = subb_u64(t.v[0], FN_N.v[0], &borrow);
        t.v[1] = subb_u64(t.v[1], FN_N.v[1], &borrow);
        t.v[2] = subb_u64(t.v[2], FN_N.v[2], &borrow);
        t.v[3] = subb_u64(t.v[3], FN_N.v[3], &borrow);
    }

    *r = t;
}

void fn_sub(fn_t *r, const fn_t *a, const fn_t *b) {
    fn_t t;
    uint64_t borrow = 0;
    fn_sub_raw(&t, a, b, &borrow);

    /* if underflow, add back n */
    if (borrow) {
        fn_add_raw(&t, &t, &FN_N);
        /* result is < 2n, reduce once */
        if (fn_cmp(&t, &FN_N) >= 0) {
            uint64_t br = 0;
            t.v[0] = subb_u64(t.v[0], FN_N.v[0], &br);
            t.v[1] = subb_u64(t.v[1], FN_N.v[1], &br);
            t.v[2] = subb_u64(t.v[2], FN_N.v[2], &br);
            t.v[3] = subb_u64(t.v[3], FN_N.v[3], &br);
        }
    }
    *r = t;
}

void fn_mul(fn_t *r, const fn_t *a, const fn_t *b) {
    /* Ensure inputs normalized for best behavior */
    fn_t aa = *a, bb = *b;
    fn_normalize(&aa);
    fn_normalize(&bb);

    u512_t x;
    mul_4x4(&x, &aa, &bb);
    barrett_mod_n(r, &x);
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

int fn_compute_t(uint8_t t_be[32],
                 const uint8_t d_be[32],
                 const uint8_t r_be[32],
                 const fn_t *x)
{
    fn_t d, rr, xr, t;

    fn_from_be(&d, d_be);
    fn_from_be(&rr, r_be);

    fn_normalize(&d);
    fn_normalize(&rr);

    if (fn_is_zero(&d) || fn_is_zero(&rr)) return 0;

    fn_mul(&xr, x, &rr);
    fn_add(&t, &d, &xr);

    if (fn_is_zero(&t)) return 0;

    fn_to_be(t_be, &t);
    return 1;
}
