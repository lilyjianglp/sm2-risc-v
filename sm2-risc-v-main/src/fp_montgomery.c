/*
 * fp_montgomery.c — 基于 Montgomery 表示的 SM2 素域 Fp 实现
 *
 * 说明：
 *   本文件仅提供 Montgomery 域内部算术与边界转换接口。
 *   若 a_bar = aR mod p，则：
 *
 *     fp_mont_add(r, a_bar, b_bar) => (a+b)R mod p
 *     fp_mont_sub(r, a_bar, b_bar) => (a-b)R mod p
 *     fp_mont_mul(r, a_bar, b_bar) => (ab)R mod p
 *     fp_mont_sqr(r, a_bar)        => (a^2)R mod p
 *     fp_mont_inv(r, a_bar)        => (a^{-1})R mod p
 *
 *   边界转换：
 *     fp_to_mont(r, a)   : a      -> aR mod p
 *     fp_from_mont(r, a) : aR mod p -> a
 *
 * 注意：
 *   本文件不实现普通域 fp_mul/fp_sqr/fp_inv/fp_reduce；
 *   普通域接口由另一份“素数专用版 fp.c”提供。
 */

#include "fp.h"
#include <stdint.h>

/* ============================================================
 *  常量
 * ============================================================ */

const fp_t FP_P = {
    .v = {
        0xFFFFFFFFFFFFFFFFULL,
        0xFFFFFFFF00000000ULL,
        0xFFFFFFFFFFFFFFFFULL,
        0xFFFFFFFEFFFFFFFFULL
    }
};

const fp_t FP_ZERO = { .v = {0, 0, 0, 0} };
const fp_t FP_ONE  = { .v = {1, 0, 0, 0} };

/* R mod p = 2^256 mod p */
const fp_t FP_MONT_ONE = {
    .v = {
        0x0000000000000001ULL,
        0x00000000FFFFFFFFULL,
        0x0000000000000000ULL,
        0x0000000100000000ULL
    }
};

/* R^2 mod p */
const fp_t FP_MONT_R2 = {
    .v = {
        0x0000000200000003ULL,
        0x00000002FFFFFFFFULL,
        0x0000000100000001ULL,
        0x0000000400000002ULL
    }
};

/* ============================================================
 *  内部工具
 * ============================================================ */

static inline uint64_t ct_mask_u64(uint64_t bit01) {
    return 0ULL - (bit01 & 1ULL);
}

static inline uint64_t addc_u64(uint64_t a, uint64_t b, uint64_t *carry) {
    uint64_t c  = *carry;
    uint64_t s  = a + b;
    uint64_t c1 = (s < a);
    uint64_t s2 = s + c;
    uint64_t c2 = (s2 < s);
    *carry = c1 | c2;
    return s2;
}

static inline uint64_t subb_u64(uint64_t a, uint64_t b, uint64_t *borrow) {
    uint64_t c  = *borrow;
    uint64_t d  = a - b;
    uint64_t b1 = (a < b);
    uint64_t d2 = d - c;
    uint64_t b2 = (d < c);
    *borrow = b1 | b2;
    return d2;
}

/* [OPT-1] 64x64 -> 128 */
static inline void mul64(uint64_t a, uint64_t b, uint64_t *lo, uint64_t *hi) {
#if defined(__SIZEOF_INT128__)
    __uint128_t z = (__uint128_t)a * b;
    *lo = (uint64_t)z;
    *hi = (uint64_t)(z >> 64);
#else
    uint64_t a0=(uint32_t)a, a1=a>>32, b0=(uint32_t)b, b1=b>>32;
    uint64_t p00=a0*b0, p01=a0*b1, p10=a1*b0, p11=a1*b1;
    uint64_t mid=(p00>>32)+(uint32_t)p01+(uint32_t)p10;
    uint64_t mid_hi=(mid>>32)+(p01>>32)+(p10>>32);
    *lo=(p00&0xFFFFFFFFULL)|(mid<<32);
    *hi=p11+mid_hi;
#endif
}

/* 192-bit accumulator: acc += (hi:lo) */
static inline void acc_add(uint64_t *a0, uint64_t *a1, uint64_t *a2,
                           uint64_t lo, uint64_t hi) {
    uint64_t c = 0;
    *a0 = addc_u64(*a0, lo, &c);
    *a1 = addc_u64(*a1, hi, &c);
    *a2 += c;
}

/* 192-bit accumulator: acc += 2*(hi:lo) */
static inline void acc_add2(uint64_t *a0, uint64_t *a1, uint64_t *a2,
                            uint64_t lo, uint64_t hi) {
    uint64_t top = hi >> 63;
    acc_add(a0, a1, a2, lo << 1, (hi << 1) | (lo >> 63));
    *a2 += top;
}

/* 无分支条件减 p */
static inline void fp_cond_sub_p(fp_t *r) {
    fp_t t;
    uint64_t borrow = 0;

    t.v[0] = subb_u64(r->v[0], FP_P.v[0], &borrow);
    t.v[1] = subb_u64(r->v[1], FP_P.v[1], &borrow);
    t.v[2] = subb_u64(r->v[2], FP_P.v[2], &borrow);
    t.v[3] = subb_u64(r->v[3], FP_P.v[3], &borrow);

    uint64_t m = ct_mask_u64(borrow ^ 1ULL);

    r->v[0] = (t.v[0] & m) | (r->v[0] & ~m);
    r->v[1] = (t.v[1] & m) | (r->v[1] & ~m);
    r->v[2] = (t.v[2] & m) | (r->v[2] & ~m);
    r->v[3] = (t.v[3] & m) | (r->v[3] & ~m);
}

/* ============================================================
 *  基础公开工具
 * ============================================================ */

int fp_cmp(const fp_t *a, const fp_t *b) {
    uint64_t gt = 0, lt = 0;
    for (int i = FP_LIMBS - 1; i >= 0; --i) {
        uint64_t u = ~(gt | lt) & 1ULL;
        gt |= (uint64_t)(a->v[i] > b->v[i]) & u;
        lt |= (uint64_t)(a->v[i] < b->v[i]) & u;
    }
    return (int)gt - (int)lt;
}

void fp_copy(fp_t *r, const fp_t *a) { *r = *a; }
void fp_set_zero(fp_t *r)            { *r = FP_ZERO; }
void fp_set_one(fp_t *r)             { *r = FP_ONE; }

int fp_is_zero(const fp_t *a) {
    return (a->v[0] | a->v[1] | a->v[2] | a->v[3]) == 0;
}

int fp_is_equal(const fp_t *a, const fp_t *b) {
    return ((a->v[0]^b->v[0])|(a->v[1]^b->v[1])|(a->v[2]^b->v[2])|(a->v[3]^b->v[3])) == 0;
}

void fp_from_bytes(fp_t *r, const uint8_t in[32]) {
    fp_t t;
    for (int i = 0; i < 4; i++) {
        uint64_t w = 0;
        for (int j = 0; j < 8; j++) w = (w << 8) | (uint64_t)in[i*8+j];
        t.v[3-i] = w;
    }
    fp_t tp; uint64_t borrow = 0;
    tp.v[0] = subb_u64(t.v[0], FP_P.v[0], &borrow);
    tp.v[1] = subb_u64(t.v[1], FP_P.v[1], &borrow);
    tp.v[2] = subb_u64(t.v[2], FP_P.v[2], &borrow);
    tp.v[3] = subb_u64(t.v[3], FP_P.v[3], &borrow);
    uint64_t m = ct_mask_u64(borrow ^ 1ULL);
    t.v[0] = (tp.v[0] & m) | (t.v[0] & ~m);
    t.v[1] = (tp.v[1] & m) | (t.v[1] & ~m);
    t.v[2] = (tp.v[2] & m) | (t.v[2] & ~m);
    t.v[3] = (tp.v[3] & m) | (t.v[3] & ~m);
    *r = t;
}

void fp_to_bytes(uint8_t out[32], const fp_t *a) {
    for (int limb = 0; limb < 4; limb++) {
        uint64_t w = a->v[3-limb];
        for (int j = 7; j >= 0; j--) out[limb*8+(7-j)] = (uint8_t)(w >> (j*8));
    }
}

/* ---- 普通域 add/sub/neg（sm2_curve.c 的 sm2_jacobian_to_affine 等用到）---- */
#ifndef USE_RV64_ASM_ADD
void fp_add(fp_t *r, const fp_t *a, const fp_t *b) {
    uint64_t carry = 0; fp_t t;
    t.v[0] = addc_u64(a->v[0], b->v[0], &carry);
    t.v[1] = addc_u64(a->v[1], b->v[1], &carry);
    t.v[2] = addc_u64(a->v[2], b->v[2], &carry);
    t.v[3] = addc_u64(a->v[3], b->v[3], &carry);
    uint64_t borrow = 0; fp_t u;
    u.v[0] = subb_u64(t.v[0], FP_P.v[0], &borrow);
    u.v[1] = subb_u64(t.v[1], FP_P.v[1], &borrow);
    u.v[2] = subb_u64(t.v[2], FP_P.v[2], &borrow);
    u.v[3] = subb_u64(t.v[3], FP_P.v[3], &borrow);
    uint64_t m = ct_mask_u64((carry & 1ULL) | ((borrow ^ 1ULL) & 1ULL));
    r->v[0]=(u.v[0]&m)|(t.v[0]&~m); r->v[1]=(u.v[1]&m)|(t.v[1]&~m);
    r->v[2]=(u.v[2]&m)|(t.v[2]&~m); r->v[3]=(u.v[3]&m)|(t.v[3]&~m);
}
#endif

#ifndef USE_RV64_ASM_SUB
void fp_sub(fp_t *r, const fp_t *a, const fp_t *b) {
    uint64_t borrow = 0; fp_t t;
    t.v[0] = subb_u64(a->v[0], b->v[0], &borrow);
    t.v[1] = subb_u64(a->v[1], b->v[1], &borrow);
    t.v[2] = subb_u64(a->v[2], b->v[2], &borrow);
    t.v[3] = subb_u64(a->v[3], b->v[3], &borrow);
    uint64_t carry = 0; fp_t u;
    u.v[0] = addc_u64(t.v[0], FP_P.v[0], &carry);
    u.v[1] = addc_u64(t.v[1], FP_P.v[1], &carry);
    u.v[2] = addc_u64(t.v[2], FP_P.v[2], &carry);
    u.v[3] = addc_u64(t.v[3], FP_P.v[3], &carry);
    uint64_t m = ct_mask_u64(borrow);
    r->v[0]=(u.v[0]&m)|(t.v[0]&~m); r->v[1]=(u.v[1]&m)|(t.v[1]&~m);
    r->v[2]=(u.v[2]&m)|(t.v[2]&~m); r->v[3]=(u.v[3]&m)|(t.v[3]&~m);
}
#endif

/* ============================================================
 *  Montgomery 域上的模加 / 模减 / 取负
 *
 *  注意：
 *    Montgomery 域中的加减法与普通域公式形式相同，因为：
 *      aR ± bR = (a ± b)R mod p
 * ============================================================ */

void fp_mont_set_zero(fp_t *r) { *r = FP_ZERO; }
void fp_mont_set_one(fp_t *r)  { *r = FP_MONT_ONE; }


#ifndef USE_RV64_ASM_ADD
void fp_mont_add(fp_t *r, const fp_t *a, const fp_t *b) {
    uint64_t carry = 0;
    fp_t t;

    t.v[0] = addc_u64(a->v[0], b->v[0], &carry);
    t.v[1] = addc_u64(a->v[1], b->v[1], &carry);
    t.v[2] = addc_u64(a->v[2], b->v[2], &carry);
    t.v[3] = addc_u64(a->v[3], b->v[3], &carry);

    uint64_t borrow = 0;
    fp_t u;

    u.v[0] = subb_u64(t.v[0], FP_P.v[0], &borrow);
    u.v[1] = subb_u64(t.v[1], FP_P.v[1], &borrow);
    u.v[2] = subb_u64(t.v[2], FP_P.v[2], &borrow);
    u.v[3] = subb_u64(t.v[3], FP_P.v[3], &borrow);

    uint64_t m = ct_mask_u64((carry & 1ULL) | ((borrow ^ 1ULL) & 1ULL));

    r->v[0] = (u.v[0] & m) | (t.v[0] & ~m);
    r->v[1] = (u.v[1] & m) | (t.v[1] & ~m);
    r->v[2] = (u.v[2] & m) | (t.v[2] & ~m);
    r->v[3] = (u.v[3] & m) | (t.v[3] & ~m);
}

#endif
#ifndef USE_RV64_ASM_SUB
void fp_mont_sub(fp_t *r, const fp_t *a, const fp_t *b) {
    uint64_t borrow = 0;
    fp_t t;

    t.v[0] = subb_u64(a->v[0], b->v[0], &borrow);
    t.v[1] = subb_u64(a->v[1], b->v[1], &borrow);
    t.v[2] = subb_u64(a->v[2], b->v[2], &borrow);
    t.v[3] = subb_u64(a->v[3], b->v[3], &borrow);

    uint64_t carry = 0;
    fp_t u;

    u.v[0] = addc_u64(t.v[0], FP_P.v[0], &carry);
    u.v[1] = addc_u64(t.v[1], FP_P.v[1], &carry);
    u.v[2] = addc_u64(t.v[2], FP_P.v[2], &carry);
    u.v[3] = addc_u64(t.v[3], FP_P.v[3], &carry);
    (void)carry;

    uint64_t m = ct_mask_u64(borrow);

    r->v[0] = (u.v[0] & m) | (t.v[0] & ~m);
    r->v[1] = (u.v[1] & m) | (t.v[1] & ~m);
    r->v[2] = (u.v[2] & m) | (t.v[2] & ~m);
    r->v[3] = (u.v[3] & m) | (t.v[3] & ~m);
}

#endif

#ifndef USE_RV64_ASM_NEG
void fp_mont_neg(fp_t *r, const fp_t *a) {
    uint64_t borrow = 0;
    fp_t t;

    t.v[0] = subb_u64(FP_P.v[0], a->v[0], &borrow);
    t.v[1] = subb_u64(FP_P.v[1], a->v[1], &borrow);
    t.v[2] = subb_u64(FP_P.v[2], a->v[2], &borrow);
    t.v[3] = subb_u64(FP_P.v[3], a->v[3], &borrow);

    uint64_t m = ct_mask_u64((uint64_t)fp_is_zero(a));

    r->v[0] = t.v[0] & ~m;
    r->v[1] = t.v[1] & ~m;
    r->v[2] = t.v[2] & ~m;
    r->v[3] = t.v[3] & ~m;
}
#endif
/* ============================================================
 *  4x4 Comba 乘法 / 平方
 * ============================================================ */

void fp_mul_4x4_u512_ct(uint64_t t[8], const fp_t *a, const fp_t *b) {
    const uint64_t a0=a->v[0], a1=a->v[1], a2=a->v[2], a3=a->v[3];
    const uint64_t b0=b->v[0], b1=b->v[1], b2=b->v[2], b3=b->v[3];
    uint64_t c0=0, c1=0, c2=0, lo, hi;

#define MA(x,y) mul64((x),(y),&lo,&hi); acc_add(&c0,&c1,&c2,lo,hi)

    c0=0;c1=0;c2=0;
    MA(a0,b0);
    t[0]=c0; c0=c1; c1=c2; c2=0;

    MA(a0,b1); MA(a1,b0);
    t[1]=c0; c0=c1; c1=c2; c2=0;

    MA(a0,b2); MA(a1,b1); MA(a2,b0);
    t[2]=c0; c0=c1; c1=c2; c2=0;

    MA(a0,b3); MA(a1,b2); MA(a2,b1); MA(a3,b0);
    t[3]=c0; c0=c1; c1=c2; c2=0;

    MA(a1,b3); MA(a2,b2); MA(a3,b1);
    t[4]=c0; c0=c1; c1=c2; c2=0;

    MA(a2,b3); MA(a3,b2);
    t[5]=c0; c0=c1; c1=c2; c2=0;

    MA(a3,b3);
    t[6]=c0; t[7]=c1;

#undef MA
}

static void comba_sqr(uint64_t t[8], const fp_t *a) {
    const uint64_t a0=a->v[0], a1=a->v[1], a2=a->v[2], a3=a->v[3];
    uint64_t c0=0, c1=0, c2=0, lo, hi;

#define MAS(x,y) mul64((x),(y),&lo,&hi); acc_add2(&c0,&c1,&c2,lo,hi)
#define MAD(x,y) mul64((x),(y),&lo,&hi); acc_add(&c0,&c1,&c2,lo,hi)

    c0=0;c1=0;c2=0;
    MAD(a0,a0);
    t[0]=c0; c0=c1; c1=c2; c2=0;

    MAS(a0,a1);
    t[1]=c0; c0=c1; c1=c2; c2=0;

    MAS(a0,a2); MAD(a1,a1);
    t[2]=c0; c0=c1; c1=c2; c2=0;

    MAS(a0,a3); MAS(a1,a2);
    t[3]=c0; c0=c1; c1=c2; c2=0;

    MAS(a1,a3); MAD(a2,a2);
    t[4]=c0; c0=c1; c1=c2; c2=0;

    MAS(a2,a3);
    t[5]=c0; c0=c1; c1=c2; c2=0;

    MAD(a3,a3);
    t[6]=c0; t[7]=c1;

#undef MAS
#undef MAD
}

/* ============================================================
 *  CIOS Montgomery reduction
 *
 *  r = t * R^{-1} mod p
 *  对 SM2 p 有 p' = 1，可取 m = t[i]
 * ============================================================ */

static void fp_mont_reduce(fp_t *r, uint64_t t[9]) {
    for (int i = 0; i < 4; ++i) {
        uint64_t m = t[i];  /* p' = 1 */
        uint64_t carry = 0;

#define CIOS_STEP(pv, ti) do {                                \
    __uint128_t _p = (__uint128_t)m * (pv) + (ti) + carry;   \
    (ti)   = (uint64_t)_p;                                   \
    carry  = (uint64_t)(_p >> 64);                           \
} while (0)

        CIOS_STEP(FP_P.v[0], t[i+0]);
        CIOS_STEP(FP_P.v[1], t[i+1]);
        CIOS_STEP(FP_P.v[2], t[i+2]);
        CIOS_STEP(FP_P.v[3], t[i+3]);

#undef CIOS_STEP

        for (int k = i + 4; k <= 8 && carry; ++k) {
            uint64_t c2 = 0;
            t[k] = addc_u64(t[k], carry, &c2);
            carry = c2;
        }
    }

    r->v[0] = t[4];
    r->v[1] = t[5];
    r->v[2] = t[6];
    r->v[3] = t[7];

    /* t[8] 产生一个额外的 2^256 倍，需要补偿 2^256 mod p = FP_MONT_ONE */
    uint64_t m_carry = ct_mask_u64(t[8] & 1ULL);
    uint64_t carry2 = 0;

    r->v[0] = addc_u64(r->v[0], FP_MONT_ONE.v[0] & m_carry, &carry2);
    r->v[1] = addc_u64(r->v[1], FP_MONT_ONE.v[1] & m_carry, &carry2);
    r->v[2] = addc_u64(r->v[2], FP_MONT_ONE.v[2] & m_carry, &carry2);
    r->v[3] = addc_u64(r->v[3], FP_MONT_ONE.v[3] & m_carry, &carry2);
    (void)carry2;

    /* 保守做法：有限次条件减 */
    fp_cond_sub_p(r);
    fp_cond_sub_p(r);
}

/* 导出给 benchmark / 测试 */
void fp_mont_reduce_ct(fp_t *r, uint64_t t[9]) {
    fp_mont_reduce(r, t);
}

/* ============================================================
 *  Montgomery 域乘法 / 平方
 * ============================================================ */

/*
 * REDC 型乘法原语：
 *   r = REDC(a * b) = a*b*R^{-1} mod p
 *
 * 只有当 a,b 已在 Montgomery 域时，
 * 它才对应通常意义上的 Montgomery 域乘法。
 */
void fp_mont_mul_ct(fp_t *r, const fp_t *a, const fp_t *b) {
    uint64_t prod[8], t[9];
    fp_mul_4x4_u512_ct(prod, a, b);
    for (int i = 0; i < 8; ++i) t[i] = prod[i];
    t[8] = 0;
    fp_mont_reduce(r, t);
}

/* 真正给上层用的 Montgomery 域乘法 */



#ifndef USE_RV64_ASM_MUL
void fp_mont_mul(fp_t *r, const fp_t *a, const fp_t *b) {
    fp_mont_mul_ct(r, a, b);
}
#endif




#ifndef USE_RV64_ASM_SQR
void fp_mont_sqr(fp_t *r, const fp_t *a) {
    uint64_t prod[8], t[9];
    comba_sqr(prod, a);
    for (int i = 0; i < 8; ++i) t[i] = prod[i];
    t[8] = 0;
    fp_mont_reduce(r, t);
}
#endif



/* ============================================================
 *  Montgomery / 普通域 边界转换
 * ============================================================ */
#ifndef USE_RV64_ASM_TO_FROM
void fp_to_mont(fp_t *r, const fp_t *a) {
    /* a -> aR mod p = REDC(a * R^2) */
    fp_mont_mul(r, a, &FP_MONT_R2);
}

void fp_from_mont(fp_t *r, const fp_t *a) {
    /* aR -> a = REDC(aR * 1) */
    fp_mont_mul(r, a, &FP_ONE);
}
#endif
/* ============================================================
 *  Montgomery 域逆元
 *
 *  若 a_bar = aR mod p，则输出 a^{-1}R mod p
 *  因为乘法/平方始终保持 Mont 域，故加法链可直接复用。
 * ============================================================ */

#define MONT_SQRN(dst, src, n) do {     \
    fp_mont_sqr((dst), (src));          \
    for (int _i = 1; _i < (n); ++_i)    \
        fp_mont_sqr((dst), (dst));      \
} while (0)

void fp_mont_inv(fp_t *r, const fp_t *a) {
    fp_t b2, b4, b8, b16, b32;
    fp_t b6, b14, b30, b31, b64;
    fp_t acc, tmp;

    /* 预计算 */
    fp_mont_sqr(&b2, a);
    fp_mont_mul(&b2, &b2, a);              /* b2 = a^(2^2-1) */

    MONT_SQRN(&b4, &b2, 2);
    fp_mont_mul(&b4, &b4, &b2);            /* b4 = a^(2^4-1) */

    MONT_SQRN(&b8, &b4, 4);
    fp_mont_mul(&b8, &b8, &b4);            /* b8 = a^(2^8-1) */

    MONT_SQRN(&b16, &b8, 8);
    fp_mont_mul(&b16, &b16, &b8);          /* b16 = a^(2^16-1) */

    MONT_SQRN(&b32, &b16, 16);
    fp_mont_mul(&b32, &b32, &b16);         /* b32 = a^(2^32-1) */

    MONT_SQRN(&b6, &b4, 2);
    fp_mont_mul(&b6, &b6, &b2);            /* b6 = a^(2^6-1) */

    MONT_SQRN(&b14, &b8, 6);
    fp_mont_mul(&b14, &b14, &b6);          /* b14 = a^(2^14-1) */

    MONT_SQRN(&b30, &b16, 14);
    fp_mont_mul(&b30, &b30, &b14);         /* b30 = a^(2^30-1) */

    fp_mont_sqr(&b31, &b30);
    fp_mont_mul(&b31, &b31, a);            /* b31 = a^(2^31-1) */

    MONT_SQRN(&b64, &b32, 32);
    fp_mont_mul(&b64, &b64, &b32);         /* b64 = a^(2^64-1) */

    /* 从 p-2 的 64-bit 分块结构组装指数 */

    /* v[3] = 0xFFFFFFFEFFFFFFFF = b32 + b31*2^33 */
    MONT_SQRN(&tmp, &b31, 33);
    fp_mont_mul(&acc, &b32, &tmp);

    /* v[2] = 0xFFFFFFFFFFFFFFFF = b64 */
    MONT_SQRN(&acc, &acc, 64);
    fp_mont_mul(&acc, &acc, &b64);

    /* v[1] = 0xFFFFFFFF00000000 */
    MONT_SQRN(&acc, &acc, 32);
    fp_mont_mul(&acc, &acc, &b32);
    MONT_SQRN(&acc, &acc, 32);

    /* v[0] = 0xFFFFFFFFFFFFFFFD = b32*2^32 + b30*2^2 + 1 */
    MONT_SQRN(&acc, &acc, 32);
    fp_mont_mul(&acc, &acc, &b32);

    MONT_SQRN(&acc, &acc, 30);
    fp_mont_mul(&acc, &acc, &b30);

    fp_mont_sqr(&acc, &acc);
    fp_mont_sqr(&acc, &acc);
    fp_mont_mul(r, &acc, a);
}

#undef MONT_SQRN