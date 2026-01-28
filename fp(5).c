#include "fp.h"
#include <string.h>
#include <stdint.h>


/* ============================================================
 *  内部工具：64-bit 加/减 带进位/借位
 * ============================================================ */

static inline uint64_t addc_u64(uint64_t a, uint64_t b, uint64_t* carry) {
    uint64_t c = *carry;          /* carry-in: 0/1 */
    uint64_t sum = a + b;
    uint64_t carry1 = (sum < a);  /* a+b 溢出则进位 */
    uint64_t sum2 = sum + c;
    uint64_t carry2 = (sum2 < sum); /* 再加 carry-in 溢出则进位 */
    *carry = carry1 | carry2;
    return sum2;
}

static inline uint64_t subb_u64(uint64_t a, uint64_t b, uint64_t* borrow) {
    uint64_t c = *borrow;          /* borrow-in: 0/1 */
    uint64_t diff = a - b;
    uint64_t borrow1 = (a < b);    /* a-b 不够减则借位 */
    uint64_t diff2 = diff - c;
    uint64_t borrow2 = (diff < c); /* 再减 borrow-in 不够减则借位 */
    *borrow = borrow1 | borrow2;
    return diff2;
}

/* ============================================================
 *  SM2 素数 p
 *  p = FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF
 *
 *  limb(小端)：
 *    v0 = 0xFFFFFFFFFFFFFFFF
 *    v1 = 0xFFFFFFFF00000000
 *    v2 = 0xFFFFFFFFFFFFFFFF
 *    v3 = 0xFFFFFFFEFFFFFFFF
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

/* ============================================================
 *  256-bit 比较：a ? b
 * ============================================================ */
int fp_cmp(const fp_t* a, const fp_t* b) {
    /* gt/lt 最终会是 0 或 1：
     * gt=1 表示 a>b；lt=1 表示 a<b；都为 0 表示相等
     */
    uint64_t gt = 0;
    uint64_t lt = 0;

    /* 从最高 limb 到最低 limb 扫描 */
    for (int i = FP_LIMBS - 1; i >= 0; --i) {
        uint64_t ai = a->v[i];
        uint64_t bi = b->v[i];

        /* 这两个比较通常会被编译成无分支 setcc/sltu */
        uint64_t ai_gt_bi = (ai > bi);  /* 0/1 */
        uint64_t ai_lt_bi = (ai < bi);  /* 0/1 */

        /* 只在“还没分出大小”(gt|lt==0) 的情况下，记录第一次出现的大小关系 */
        uint64_t undecided = ~(gt | lt) & 1ULL;
        gt |= ai_gt_bi & undecided;
        lt |= ai_lt_bi & undecided;
    }

    /* 把 0/1 映射到 -1/0/1 */
    return (int)gt - (int)lt;
}

void fp_copy(fp_t* r, const fp_t* a) {
    *r = *a;
}

void fp_set_zero(fp_t* r) {
    r->v[0] = r->v[1] = r->v[2] = r->v[3] = 0;
}

void fp_set_one(fp_t* r) {
    r->v[0] = 1;
    r->v[1] = r->v[2] = r->v[3] = 0;
}

int fp_is_zero(const fp_t* a) {
    return (a->v[0] | a->v[1] | a->v[2] | a->v[3]) == 0;
}

int fp_is_equal(const fp_t* a, const fp_t* b) {
    uint64_t x = 0;
    x |= (a->v[0] ^ b->v[0]);
    x |= (a->v[1] ^ b->v[1]);
    x |= (a->v[2] ^ b->v[2]);
    x |= (a->v[3] ^ b->v[3]);
    return x == 0;
}


/* ============================================================
 *  模加 / 模减 / 取负
 * ============================================================ */

static inline uint64_t ct_mask_u64(uint64_t bit01) {
    /* bit01 must be 0 or 1 -> returns 0x00..00 or 0xFF..FF */
    return 0ULL - (bit01 & 1ULL);
}



#ifndef USE_ASM_FP_ADD
void fp_add(fp_t* r, const fp_t* a, const fp_t* b) {
    uint64_t carry = 0;
    fp_t t;

    t.v[0] = addc_u64(a->v[0], b->v[0], &carry);
    t.v[1] = addc_u64(a->v[1], b->v[1], &carry);
    t.v[2] = addc_u64(a->v[2], b->v[2], &carry);
    t.v[3] = addc_u64(a->v[3], b->v[3], &carry);

    /* Always compute u = t - p */
    uint64_t borrow = 0;
    fp_t u;
    u.v[0] = subb_u64(t.v[0], FP_P.v[0], &borrow);
    u.v[1] = subb_u64(t.v[1], FP_P.v[1], &borrow);
    u.v[2] = subb_u64(t.v[2], FP_P.v[2], &borrow);
    u.v[3] = subb_u64(t.v[3], FP_P.v[3], &borrow);

    /* do_sub = carry || (t >= p)
       t >= p <=> borrow == 0
       so do_sub = carry | (1 ^ borrow)
    */
    uint64_t do_sub = (carry & 1ULL) | ((borrow ^ 1ULL) & 1ULL);
    uint64_t m = ct_mask_u64(do_sub);

    /* r = do_sub ? u : t  (branchless) */
    r->v[0] = (u.v[0] & m) | (t.v[0] & ~m);
    r->v[1] = (u.v[1] & m) | (t.v[1] & ~m);
    r->v[2] = (u.v[2] & m) | (t.v[2] & ~m);
    r->v[3] = (u.v[3] & m) | (t.v[3] & ~m);
}
#endif


#ifndef USE_ASM_FP_SUB
void fp_sub(fp_t* r, const fp_t* a, const fp_t* b) {
    uint64_t borrow = 0;
    fp_t t;

    t.v[0] = subb_u64(a->v[0], b->v[0], &borrow);
    t.v[1] = subb_u64(a->v[1], b->v[1], &borrow);
    t.v[2] = subb_u64(a->v[2], b->v[2], &borrow);
    t.v[3] = subb_u64(a->v[3], b->v[3], &borrow);

    /* Always compute u = t + p */
    uint64_t carry = 0;
    fp_t u;
    u.v[0] = addc_u64(t.v[0], FP_P.v[0], &carry);
    u.v[1] = addc_u64(t.v[1], FP_P.v[1], &carry);
    u.v[2] = addc_u64(t.v[2], FP_P.v[2], &carry);
    u.v[3] = addc_u64(t.v[3], FP_P.v[3], &carry);
    (void)carry; /* carry 可忽略：t < p 时 t+p < 2p < 2^256（对常见 256-bit 素数模数成立） */

    /* r = borrow ? u : t (branchless) */
    uint64_t m = ct_mask_u64(borrow);

    r->v[0] = (u.v[0] & m) | (t.v[0] & ~m);
    r->v[1] = (u.v[1] & m) | (t.v[1] & ~m);
    r->v[2] = (u.v[2] & m) | (t.v[2] & ~m);
    r->v[3] = (u.v[3] & m) | (t.v[3] & ~m);
}

#endif


#ifndef USE_ASM_FP_NEG
void fp_neg(fp_t* r, const fp_t* a) {
    /* t = p - a (always) */
    uint64_t borrow = 0;
    fp_t t;
    t.v[0] = subb_u64(FP_P.v[0], a->v[0], &borrow);
    t.v[1] = subb_u64(FP_P.v[1], a->v[1], &borrow);
    t.v[2] = subb_u64(FP_P.v[2], a->v[2], &borrow);
    t.v[3] = subb_u64(FP_P.v[3], a->v[3], &borrow);

    /* if a == 0 then r = 0 else r = t */
    uint64_t is_zero = fp_is_zero(a);     // 1 if a==0 else 0
    uint64_t m = ct_mask_u64(is_zero);       // all-ones if a==0

    r->v[0] = t.v[0] & ~m;
    r->v[1] = t.v[1] & ~m;
    r->v[2] = t.v[2] & ~m;
    r->v[3] = t.v[3] & ~m;
}
#endif


/* ============================================================
 *  512-bit 工具（用于模约减）
 * ============================================================ */

/* ============================================================
 *  SM2 专用素数约减（路线B）：512-bit -> 256-bit
 *
 *  p = 2^256 - 2^224 - 2^96 + 2^64 - 1
 *  => 2^256 ≡ 1 + 2^224 + 2^96 - 2^64 (mod p)
 *
 *  直接把 t4..t7 的贡献展开折回低位，避免 many-round fold。
 * ============================================================ */

/* ============================================================
 *  SM2 专用素数约减（路线B）：512-bit -> 256-bit
 *  去除 __int128：用 (lo, hi_units) 显式 carry/borrow 传播
 * ============================================================ */

typedef struct {
    uint64_t lo[6];   /* limb 本体（mod 2^64） */
    int64_t  hi[6];   /* 以 2^64 为单位的“有符号进位/借位”，最终要传播到更高 limb */
} acc6_t;

static inline void acc6_zero(acc6_t *a) {
    for (int i = 0; i < 6; i++) { a->lo[i] = 0; a->hi[i] = 0; }
}

/* lo[idx] += v; hi[idx] += carry(0/1) */
static inline void acc6_add_u64(acc6_t *a, int idx, uint64_t v) {
    uint64_t old = a->lo[idx];
    uint64_t sum = old + v;
    a->lo[idx] = sum;
    a->hi[idx] += (int64_t)(sum < old);
}

/* lo[idx] -= v; hi[idx] -= borrow(0/1) */
static inline void acc6_sub_u64(acc6_t *a, int idx, uint64_t v) {
    uint64_t old = a->lo[idx];
    uint64_t diff = old - v;
    a->lo[idx] = diff;
    a->hi[idx] -= (int64_t)(old < v);
}

/* constant-time: 将 signed q 加到 limb[idx]（q 以“普通整数”形式加到 lo 上）
 * 若 q<0 则等价于减 abs(q)。更新 hi[idx] 记录 carry/borrow。
 * 说明：这里 q 的大小通常很小（由多次 carry/borrow 累积），但实现对任意 int64 都正确。
 */
static inline void acc6_add_i64_ct(acc6_t *a, int idx, int64_t q) {
    uint64_t old = a->lo[idx];
    uint64_t uq  = (uint64_t)q;

    /* sign = 1 if q<0 else 0 */
    uint64_t sign = uq >> 63;
    uint64_t m    = 0ULL - sign;                 /* all-ones if negative */

    /* mag = abs(q) (two's complement, CT) */
    uint64_t mag  = (uq ^ m) - m;

    uint64_t sum  = old + mag;
    uint64_t diff = old - mag;

    uint64_t newlo = (diff & m) | (sum & ~m);
    uint64_t carry_add  = (sum  < old);          /* valid for add path */
    uint64_t borrow_sub = (old  < mag);          /* valid for sub path */

    a->lo[idx] = newlo;

    /* hi += carry when q>=0 ; hi -= borrow when q<0 (CT) */
    a->hi[idx] += (int64_t)(carry_add  & (sign ^ 1ULL));
    a->hi[idx] -= (int64_t)(borrow_sub & sign);
}

/* shift-add: acc += v<<sh （只用 64-bit add + carry 记到 hi） */
static inline void acc_add_shift_u64(acc6_t *a, uint64_t v, int sh) {
    int w = sh >> 6;
    int b = sh & 63;
    if ((unsigned)w >= 6) return;

    if (b == 0) {
        acc6_add_u64(a, w, v);
    } else {
        uint64_t lo = v << b;
        uint64_t hi = v >> (64 - b);
        acc6_add_u64(a, w, lo);
        if (w + 1 < 6) acc6_add_u64(a, w + 1, hi);
    }
}

/* shift-sub: acc -= v<<sh */
static inline void acc_sub_shift_u64(acc6_t *a, uint64_t v, int sh) {
    int w = sh >> 6;
    int b = sh & 63;
    if ((unsigned)w >= 6) return;

    if (b == 0) {
        acc6_sub_u64(a, w, v);
    } else {
        uint64_t lo = v << b;
        uint64_t hi = v >> (64 - b);
        acc6_sub_u64(a, w, lo);
        if (w + 1 < 6) acc6_sub_u64(a, w + 1, hi);
    }
}

/* 显式归一化：out[i]=lo[i]；把 hi[i]（以 2^64 为单位的有符号进位）传播到 i+1 的 lo 里
 * 完全不需要 __int128 / 除法；并且传播使用 CT 的 acc6_add_i64_ct。
 */
static inline void normalize_acc6(acc6_t *a, uint64_t out[6]) {
    for (int i = 0; i < 6; i++) {
        out[i] = a->lo[i];
        int64_t q = a->hi[i];

        /* 清空当前 limb（可选，但更像原先语义，避免误用） */
        a->lo[i] = 0;
        a->hi[i] = 0;

        /* 传播到下一 limb：q * 2^64 * 2^(64*i) == q * 2^(64*(i+1)) */
        if (i + 1 < 6) {
            acc6_add_i64_ct(a, i + 1, q);
        }
    }
}

/* k*2^256 折回：k*(1 + 2^224 + 2^96 - 2^64) */
static inline void fold_k_2p256(acc6_t *a, uint64_t k) {
    acc_add_shift_u64(a, k,   0);
    acc_add_shift_u64(a, k,  96);
    acc_add_shift_u64(a, k, 224);
    acc_sub_shift_u64(a, k,  64);
}

/* 2^320 ≡ 1 + 2^32 + 2^160 + 2^224 (mod p) */
static inline void fold_k_2p320(acc6_t *a, uint64_t k) {
    acc_add_shift_u64(a, k,   0);
    acc_add_shift_u64(a, k,  32);
    acc_add_shift_u64(a, k, 160);
    acc_add_shift_u64(a, k, 224);
}

/* 2^384 ≡ 1 + 2^32 + 2^96 + 2^128 + 2^225 (mod p) */
static inline void fold_k_2p384(acc6_t *a, uint64_t k) {
    acc_add_shift_u64(a, k,   0);
    acc_add_shift_u64(a, k,  32);
    acc_add_shift_u64(a, k,  96);
    acc_add_shift_u64(a, k, 128);
    acc_add_shift_u64(a, k, 225);
}

/* 2^448 ≡ 2 + 2^33 + 2^96 + 2^129 + 2^160 + 2^192 + 2^225 - 2^64 (mod p) */
static inline void fold_k_2p448(acc6_t *a, uint64_t k) {
    acc_add_shift_u64(a, k,   1);  /* 2*k */
    acc_add_shift_u64(a, k,  33);
    acc_add_shift_u64(a, k,  96);
    acc_add_shift_u64(a, k, 129);
    acc_add_shift_u64(a, k, 160);
    acc_add_shift_u64(a, k, 192);
    acc_add_shift_u64(a, k, 225);
    acc_sub_shift_u64(a, k,  64);
}


static inline void fp_cond_sub_p(fp_t* r) {
    uint64_t borrow = 0;
    fp_t t;
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

/* ✅ 512-bit T (8 limbs, little-endian) -> r (4 limbs): r = T mod p */
static void mod_reduce_512_to_256(fp_t* r, const uint64_t T[8]) {
    uint64_t t0=T[0], t1=T[1], t2=T[2], t3=T[3];
    uint64_t t4=T[4], t5=T[5], t6=T[6], t7=T[7];

    acc6_t acc;
    acc6_zero(&acc);
    acc.lo[0] = t0;
    acc.lo[1] = t1;
    acc.lo[2] = t2;
    acc.lo[3] = t3;

    fold_k_2p256(&acc, t4);
    fold_k_2p320(&acc, t5);
    fold_k_2p384(&acc, t6);
    fold_k_2p448(&acc, t7);

    uint64_t x[6];
    normalize_acc6(&acc, x);

    /* 固定 2 轮尾巴折回（避免 while，固定时间） */
    for (int round = 0; round < 2; round++) {
        uint64_t k4 = x[4];
        uint64_t k5 = x[5];
        x[4] = 0;
        x[5] = 0;

        acc6_t acc2;
        acc6_zero(&acc2);
        acc2.lo[0] = x[0];
        acc2.lo[1] = x[1];
        acc2.lo[2] = x[2];
        acc2.lo[3] = x[3];

        fold_k_2p256(&acc2, k4);
        fold_k_2p320(&acc2, k5);

        normalize_acc6(&acc2, x);
    }

    fp_t out;
    out.v[0] = x[0];
    out.v[1] = x[1];
    out.v[2] = x[2];
    out.v[3] = x[3];

    /* 更保守：固定做 3 次条件减 p（仍然 CT） */
    fp_cond_sub_p(&out);
    fp_cond_sub_p(&out);
    fp_cond_sub_p(&out);

    *r = out;
}


/* ============================================================
 *  模乘 / 平方
 * ============================================================ */
static inline void mul64x64_128(uint64_t a, uint64_t b,
                                uint64_t *lo, uint64_t *hi)
{
    uint64_t a0 = (uint32_t)a;
    uint64_t a1 = a >> 32;
    uint64_t b0 = (uint32_t)b;
    uint64_t b1 = b >> 32;

    uint64_t p00 = a0 * b0;          // 64
    uint64_t p01 = a0 * b1;          // 64
    uint64_t p10 = a1 * b0;          // 64
    uint64_t p11 = a1 * b1;          // 64

    /* lo = p00 + ((p01 + p10) << 32) */
    uint64_t mid = p01 + p10;
    uint64_t carry_mid = (mid < p01);                // mid 溢出

    uint64_t lo_ = p00 + (mid << 32);
    uint64_t carry_lo = (lo_ < p00);                 // lo 溢出

    uint64_t hi_ = p11 + (mid >> 32) + (carry_mid << 32) + carry_lo;

    *lo = lo_;
    *hi = hi_;
}


/* acc = acc + (lo + 2^64 * hi), acc 是 192-bit */
static inline void acc192_add128(uint64_t *acc0,
                                 uint64_t *acc1,
                                 uint64_t *acc2,
                                 uint64_t lo,
                                 uint64_t hi)
{
    uint64_t carry = 0;

    /* acc0 += lo */
    *acc0 = addc_u64(*acc0, lo, &carry);

    /* acc1 += hi + carry */
    *acc1 = addc_u64(*acc1, hi, &carry);

    /* acc2 += carry */
    *acc2 += carry;
}


void fp_mul(fp_t *r, const fp_t *a, const fp_t *b) {
    const uint64_t a0 = a->v[0], a1 = a->v[1], a2 = a->v[2], a3 = a->v[3];
    const uint64_t b0 = b->v[0], b1 = b->v[1], b2 = b->v[2], b3 = b->v[3];

    uint64_t t[8];

    /* carry = c0 + 2^64*c1 */
    uint64_t c0 = 0, c1 = 0;

    uint64_t acc0, acc1, acc2;
    uint64_t lo, hi;

    /* k = 0: a0*b0 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a0, b0, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    t[0] = acc0; c0 = acc1; c1 = acc2;

    /* k = 1: a0*b1 + a1*b0 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a0, b1, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a1, b0, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    t[1] = acc0; c0 = acc1; c1 = acc2;

    /* k = 2: a0*b2 + a1*b1 + a2*b0 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a0, b2, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a1, b1, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a2, b0, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    t[2] = acc0; c0 = acc1; c1 = acc2;

    /* k = 3: a0*b3 + a1*b2 + a2*b1 + a3*b0 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a0, b3, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a1, b2, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a2, b1, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a3, b0, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    t[3] = acc0; c0 = acc1; c1 = acc2;

    /* k = 4: a1*b3 + a2*b2 + a3*b1 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a1, b3, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a2, b2, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a3, b1, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    t[4] = acc0; c0 = acc1; c1 = acc2;

    /* k = 5: a2*b3 + a3*b2 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a2, b3, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a3, b2, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    t[5] = acc0; c0 = acc1; c1 = acc2;

    /* k = 6: a3*b3 */
    acc0 = c0; acc1 = c1; acc2 = 0; 
    mul64x64_128(a3, b3, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    t[6] = acc0; c0 = acc1; c1 = acc2;

    /* k = 7: only carry remains */
    t[7] = c0;
    /* c1 should be 0 for 256x256 -> 512, but even if not, we'd be beyond 512 bits */

    /* SM2 素数专用约减（你已有：固定轮数 + CT cond-sub） */
    mod_reduce_512_to_256(r, t);
}


static inline void acc192_add128_dbl(uint64_t *acc0, uint64_t *acc1, uint64_t *acc2,
                                     uint64_t lo, uint64_t hi) {
    /* 2*(hi:lo) = ( (hi<<1 | lo>>63) : (lo<<1) ) + carry128 * 2^128 */
    uint64_t carry128 = hi >> 63;
    uint64_t dlo = lo << 1;
    uint64_t dhi = (hi << 1) | (lo >> 63);
    acc192_add128(acc0, acc1, acc2, dlo, dhi);
    *acc2 += carry128;
}


void fp_sqr(fp_t* r, const fp_t* a) {
    const uint64_t a0=a->v[0], a1=a->v[1], a2=a->v[2], a3=a->v[3];

    uint64_t t[8];

    /* carry = c0 + 2^64*c1 */
    uint64_t c0 = 0, c1 = 0;

    uint64_t acc0, acc1, acc2;
    uint64_t lo, hi;

    /* k=0: a0^2 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a0, a0, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    t[0] = acc0; c0 = acc1; c1 = acc2;

    /* k=1: 2*a0*a1 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a0, a1, &lo, &hi); acc192_add128_dbl(&acc0, &acc1, &acc2, lo, hi);
    t[1] = acc0; c0 = acc1; c1 = acc2;

    /* k=2: 2*a0*a2 + a1^2 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a0, a2, &lo, &hi); acc192_add128_dbl(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a1, a1, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    t[2] = acc0; c0 = acc1; c1 = acc2;

    /* k=3: 2*a0*a3 + 2*a1*a2 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a0, a3, &lo, &hi); acc192_add128_dbl(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a1, a2, &lo, &hi); acc192_add128_dbl(&acc0, &acc1, &acc2, lo, hi);
    t[3] = acc0; c0 = acc1; c1 = acc2;

    /* k=4: 2*a1*a3 + a2^2 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a1, a3, &lo, &hi); acc192_add128_dbl(&acc0, &acc1, &acc2, lo, hi);
    mul64x64_128(a2, a2, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    t[4] = acc0; c0 = acc1; c1 = acc2;

    /* k=5: 2*a2*a3 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a2, a3, &lo, &hi); acc192_add128_dbl(&acc0, &acc1, &acc2, lo, hi);
    t[5] = acc0; c0 = acc1; c1 = acc2;

    /* k=6: a3^2 */
    acc0 = c0; acc1 = c1; acc2 = 0;
    mul64x64_128(a3, a3, &lo, &hi); acc192_add128(&acc0, &acc1, &acc2, lo, hi);
    t[6] = acc0; c0 = acc1; c1 = acc2;

    /* k=7: only carry remains */
    t[7] = c0;

    mod_reduce_512_to_256(r, t);
}

/* ============================================================
 *  逆元：Fermat a^(p-2) mod p
 *  p-2 = FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFD
 * ============================================================ */

static uint32_t get_bits_p_minus_2(int bit_lo, int width) {
    /* exponent = p - 2 (256-bit), little-endian 64-bit limbs */
    const uint64_t e[4] = {
        0xFFFFFFFFFFFFFFFDULL,
        0xFFFFFFFF00000000ULL,
        0xFFFFFFFFFFFFFFFFULL,
        0xFFFFFFFEFFFFFFFFULL
    };

    int limb = bit_lo / 64;
    int off  = bit_lo % 64;

    if (off + width <= 64) {
        return (uint32_t)((e[limb] >> off) & ((1ULL << width) - 1ULL));
    } else {
        /* 需要跨 limb 拼接（width=4 时很少跨，但写全更稳） */
        int hi_bits = off + width - 64;
        uint64_t lo = e[limb] >> off;
        uint64_t hi = (limb + 1 < 4) ? (e[limb + 1] & ((1ULL << hi_bits) - 1ULL)) : 0;
        return (uint32_t)(lo | (hi << (64 - off)));
    }
}

/* constant-time eq for 32-bit: return 1 if a==b else 0 (no branches) */
static inline uint32_t ct_eq_u32(uint32_t a, uint32_t b) {
    uint32_t x = a ^ b;           /* 0 if equal */
    x |= x >> 16;
    x |= x >> 8;
    x |= x >> 4;
    x |= x >> 2;
    x |= x >> 1;
    return (x ^ 1u) & 1u;
}

/* constant-time select: out = table[w], scanning all 16 entries */
static inline void fp_table16_select(fp_t *out, const fp_t table[16], uint32_t w /*0..15*/) {
    fp_t r;
    r.v[0] = r.v[1] = r.v[2] = r.v[3] = 0;

    for (uint32_t i = 0; i < 16; i++) {
        uint64_t m = 0 - (uint64_t)ct_eq_u32(i, w);  /* 0xFFFF.. if i==w else 0 */
        r.v[0] |= table[i].v[0] & m;
        r.v[1] |= table[i].v[1] & m;
        r.v[2] |= table[i].v[2] & m;
        r.v[3] |= table[i].v[3] & m;
    }
    *out = r;
}

void fp_inv(fp_t* r, const fp_t* a) {
    /* 调用者保证 a != 0 */

    /* 预计算表：pow[i] = a^i, i=0..15 */
    fp_t pow[16];
    fp_set_one(&pow[0]);
    fp_copy(&pow[1], a);

    fp_sqr(&pow[2], a);
    for (int i = 3; i < 16; i++) {
        fp_mul(&pow[i], &pow[i - 1], a);
    }

    fp_t result;
    fp_set_one(&result);

    /* 固定窗口 w=4：从最高 nibble 到最低 nibble */
    for (int bit_lo = 252; bit_lo >= 0; bit_lo -= 4) {
        /* result = result^16 (4 squares) */
        fp_sqr(&result, &result);
        fp_sqr(&result, &result);
        fp_sqr(&result, &result);
        fp_sqr(&result, &result);

        /* 公开常量 nibble w in [0..15] */
        uint32_t w = get_bits_p_minus_2(bit_lo, 4);

        /* 严格 CT：扫描表选出 pow[w]，每轮固定做一次乘法（包括 w=0 乘以 1） */
        fp_t mulv;
        fp_table16_select(&mulv, pow, w);
        fp_mul(&result, &result, &mulv);
    }

    fp_copy(r, &result);
}

/* ============================================================
 *  字节序列化：32-byte big-endian
 * ============================================================ */

void fp_from_bytes(fp_t* r, const uint8_t in[32]) {
    fp_t t;
    for (int i = 0; i < 4; ++i) t.v[i] = 0;

    for (int limb = 0; limb < 4; ++limb) {
        uint64_t w = 0;
        for (int j = 0; j < 8; ++j) {
            w = (w << 8) | (uint64_t)in[limb * 8 + j];
        }
        t.v[3 - limb] = w;
    }

    /* 计算 t - p（无条件），再用 borrow 选择 */
    fp_t tp;
    uint64_t borrow = 0;
    tp.v[0] = subb_u64(t.v[0], FP_P.v[0], &borrow);
    tp.v[1] = subb_u64(t.v[1], FP_P.v[1], &borrow);
    tp.v[2] = subb_u64(t.v[2], FP_P.v[2], &borrow);
    tp.v[3] = subb_u64(t.v[3], FP_P.v[3], &borrow);

    /* borrow==0 => t>=p => 选 tp；borrow==1 => t<p => 选 t */
    uint64_t m = ct_mask_u64(borrow ^ 1ULL);  /* borrow=0 -> all1, borrow=1 -> 0 */
    t.v[0] = (tp.v[0] & m) | (t.v[0] & ~m);
    t.v[1] = (tp.v[1] & m) | (t.v[1] & ~m);
    t.v[2] = (tp.v[2] & m) | (t.v[2] & ~m);
    t.v[3] = (tp.v[3] & m) | (t.v[3] & ~m);

    *r = t;
}


void fp_to_bytes(uint8_t out[32], const fp_t* a) {
    /* limbs little-endian -> big-endian bytes */
    for (int limb = 0; limb < 4; ++limb) {
        uint64_t w = a->v[3 - limb]; /* v[3] 是最高 64 位 */
        for (int j = 7; j >= 0; --j) {
            out[limb * 8 + (7 - j)] = (uint8_t)((w >> (j * 8)) & 0xFF);
        }
    }
}


