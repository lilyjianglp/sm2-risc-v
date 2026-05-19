/* ============================================================
 * sm2_scalar.c 
 * ============================================================ */
#include "sm2_scalar.h"
#include "sm2_curve.h"
#include <string.h>

/* ============================================================
 * 固定基点预计算表
 *
 * 方法：把 256-bit scalar 拆成 32 个 8-bit 窗口
 *   T_base[i][j] = (j+1) * 2^(8*i) * G,  i=0..31, j=0..254
 *   点乘 = 32 次查表 + 32 次 mixed add，0 次 double
 *
 * 内存：32 * 255 * sizeof(sm2_affine_t) ≈ 32*255*64 = 522KB
 * ============================================================ */
#define FIXEDBASE_WINDOWS  32
#define FIXEDBASE_TABLE_SZ 255   /* j=0..254，对应 1P..255P */

static sm2_affine_t T_base[FIXEDBASE_WINDOWS][FIXEDBASE_TABLE_SZ];
static int          T_base_ready = 0;

/* 初始化固定基点表（只算一次） */
static void sm2_init_fixedbase_table(void)
{
    if (T_base_ready) return;

    sm2_affine_t G;
    sm2_get_base_affine(&G);

    /* 计算 T_base[0][0..254] = 1G .. 255G */
    sm2_jacobian_t tmp[FIXEDBASE_TABLE_SZ];
    fp_t prefix_buf[FIXEDBASE_TABLE_SZ + 1];

    sm2_jacobian_t J;
    sm2_affine_to_jacobian(&J, &G);
    tmp[0] = J;
    for (int j = 1; j < FIXEDBASE_TABLE_SZ; j++) {
        sm2_add_ja(&tmp[j], &tmp[j-1], &G);
    }
    sm2_batch_normalize(T_base[0], tmp, FIXEDBASE_TABLE_SZ, prefix_buf);
    for (int j = 0; j < FIXEDBASE_TABLE_SZ; j++) T_base[0][j].infinity = 0;

    /* 计算 T_base[i][j] = 2^(8*i) * T_base[0][j]
       即每组是上一组的 2^8 倍 */
    for (int i = 1; i < FIXEDBASE_WINDOWS; i++) {
        /* T_base[i][j] = 2^8 * T_base[i-1][j] */
        for (int j = 0; j < FIXEDBASE_TABLE_SZ; j++) {
            sm2_jacobian_t cur;
            sm2_affine_to_jacobian(&cur, &T_base[i-1][j]);
            /* double 8 次 = 乘以 2^8 */
            for (int d = 0; d < 8; d++) {
                sm2_double_jm_a_minus3(&cur, &cur);
            }
            sm2_jacobian_to_affine(&T_base[i][j], &cur);
            T_base[i][j].infinity = 0;
        }
    }

    T_base_ready = 1;
}

/* 固定基点快速标量乘：result = k * G
   k[0] = MSB, k[31] = LSB
   0 次 double，最多 32 次 mixed add */
static void sm2_scalar_mul_fixedbase(sm2_jacobian_t *R, const uint8_t k[32])
{
    sm2_init_fixedbase_table();
    sm2_jacobian_set_infinity(R);

    for (int i = 0; i < FIXEDBASE_WINDOWS; i++) {
        /* k[31-i] 是第 i 个窗口（从低位窗口开始） */
        uint8_t idx = k[31 - i];
        if (idx == 0) continue;
        /* T_base[i][idx-1] = idx * 2^(8*i) * G */
        sm2_add_ja(R, R, &T_base[i][idx - 1]);
    }
}

#ifdef USE_FP_MONT
static sm2_affine_t T_base_mont[FIXEDBASE_WINDOWS][FIXEDBASE_TABLE_SZ];
static int          T_base_mont_ready = 0;

static void sm2_init_fixedbase_table_mont(void)
{
    if (T_base_mont_ready) return;

    sm2_affine_t Gm;
    sm2_get_base_affine_mont(&Gm);

    sm2_jacobian_t tmp[FIXEDBASE_TABLE_SZ];
    fp_t prefix_buf[FIXEDBASE_TABLE_SZ + 1];

    sm2_jacobian_t J;
    sm2_affine_mont_to_jacobian_mont(&J, &Gm);
    tmp[0] = J;

    for (int j = 1; j < FIXEDBASE_TABLE_SZ; j++) {
        sm2_add_ja_mont(&tmp[j], &tmp[j - 1], &Gm);
    }

    sm2_batch_normalize_mont(T_base_mont[0], tmp, FIXEDBASE_TABLE_SZ, prefix_buf);
    for (int j = 0; j < FIXEDBASE_TABLE_SZ; j++) T_base_mont[0][j].infinity = 0;

    for (int i = 1; i < FIXEDBASE_WINDOWS; i++) {
        for (int j = 0; j < FIXEDBASE_TABLE_SZ; j++) {
            sm2_jacobian_t cur;
            sm2_affine_mont_to_jacobian_mont(&cur, &T_base_mont[i - 1][j]);
            for (int d = 0; d < 8; d++) {
                sm2_double_jm_a_minus3_mont(&cur, &cur);
            }
            sm2_jacobian_mont_to_affine_mont(&T_base_mont[i][j], &cur);
            T_base_mont[i][j].infinity = 0;
        }
    }

    T_base_mont_ready = 1;
}

static void sm2_scalar_mul_fixedbase_mont(sm2_jacobian_t *R, const uint8_t k[32])
{
    sm2_init_fixedbase_table_mont();

    int started = 0;
    sm2_jacobian_set_infinity_mont(R);

    for (int i = 0; i < FIXEDBASE_WINDOWS; i++) {
        uint8_t idx = k[31 - i];
        if (idx == 0) continue;

        if (!started) {
            sm2_affine_mont_to_jacobian_mont(R, &T_base_mont[i][idx - 1]);
            started = 1;
        } else {
            sm2_add_ja_mont(R, R, &T_base_mont[i][idx - 1]);
        }
    }

    if (!started) {
        sm2_jacobian_set_infinity_mont(R);
    }
}
#endif /* USE_FP_MONT */

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
#ifdef USE_FP_MONT
static void affine_cond_neg_mont(sm2_affine_t* p, uint32_t neg)
{
    fp_t yneg;
    fp_mont_neg(&yneg, &p->y);
    fp_cmov(&p->y, &yneg, neg);
}
#endif /* USE_FP_MONT */
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

#ifdef USE_FP_MONT
static void sm2_precompute_table_1_to_16_affine_mont(sm2_affine_t T16[16],
                                                     const sm2_affine_t* P_mont)
{
    sm2_jacobian_t tmp[16];
    sm2_jacobian_t J;

    sm2_affine_mont_to_jacobian_mont(&J, P_mont);

    tmp[0] = J;

    for (int i = 1; i < 16; i++) {
        sm2_add_ja_mont(&tmp[i], &tmp[i - 1], P_mont);
    }

    fp_t prefix_buf[16 + 1];
    sm2_batch_normalize_mont(T16, tmp, 16, prefix_buf);

    for (int i = 0; i < 16; i++) T16[i].infinity = 0;
}
#endif /* USE_FP_MONT */

/* ============================================================
 *   SchemeB: fixed w=5, fully-unrolled 5 doubles, CT select
 *    Works for BOTH odd/even k now.
 *    当 P == G（基点）时，自动走固定基点快速路径（0次double）
 * ============================================================ */
void sm2_scalar_mul_window_ct_schemeB(sm2_jacobian_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32])
{
    if (sm2_affine_is_infinity(P) || k_is_zero(k)) {
        sm2_jacobian_set_infinity(R);
        return;
    }

    /* 检查 P 是否是基点 G，若是则走固定基点快速路径 */
    {
        sm2_affine_t G;
        sm2_get_base_affine(&G);
        /* 常数时间比较：比较 x 和 y 坐标的 4 个 limb */
        uint64_t diff = 0;
        for (int i = 0; i < 4; i++) {
            diff |= (P->x.v[i] ^ G.x.v[i]);
            diff |= (P->y.v[i] ^ G.y.v[i]);
        }
        if (diff == 0) {
            sm2_scalar_mul_fixedbase(R, k);
            return;
        }
    }

    /* 任意点：原始窗口法 */
    const unsigned w = 5;
    const int n = (256 + 4) / 5; /* 52 */

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
/* ============================================================
 * Montgomery SchemeB core
 *
 * 核心原则：
 * 1. core 函数只接受 Montgomery affine 输入
 * 2. core 函数只输出 Montgomery jacobian
 * 3. 不在 core 里做 affine<->mont 的来回转换
 * 4. 包装函数再按需要做输入/输出转换
 * ============================================================ */

#ifdef USE_FP_MONT
static void sm2_scalar_mul_window_ct_schemeB_mont_core(
    sm2_jacobian_t* R,
    const sm2_affine_t* P_mont,
    const uint8_t k[32])
{
    if (sm2_affine_is_infinity(P_mont) || k_is_zero(k)) {
        sm2_jacobian_set_infinity_mont(R);
        return;
    }

    /* 检查是否为 Montgomery 域基点 Gm，若是直接走 fixed-base 快速路径 */
    {
        sm2_affine_t Gm;
        sm2_get_base_affine_mont(&Gm);

        uint64_t diff = 0;
        for (int i = 0; i < 4; i++) {
            diff |= (P_mont->x.v[i] ^ Gm.x.v[i]);
            diff |= (P_mont->y.v[i] ^ Gm.y.v[i]);
        }
        if (diff == 0) {
            sm2_scalar_mul_fixedbase_mont(R, k);
            return;
        }
    }

    const unsigned w = 5;
    const int n = (256 + 4) / 5;

    sm2_affine_t T16[16];
    sm2_precompute_table_1_to_16_affine_mont(T16, P_mont);

    /* 注意：这里不再强行覆盖 T16[0]
       sm2_precompute_table_1_to_16_affine_mont 已经生成了正确的 1P..16P 表 */

    uint64_t a[4];
    scalar_be_to_u64le(a, k);

    sm2_jacobian_set_infinity_mont(R);
    int started = 0;

    for (int i = n - 1; i >= 0; i--) {
        if (started) {
            sm2_double_jm_a_minus3_mont(R, R);
            sm2_double_jm_a_minus3_mont(R, R);
            sm2_double_jm_a_minus3_mont(R, R);
            sm2_double_jm_a_minus3_mont(R, R);
            sm2_double_jm_a_minus3_mont(R, R);
        }

        int di = booth_get_digit_u64le(a, w, i);

        int32_t d32 = (int32_t)di;
        uint32_t sign = (uint32_t)(d32 >> 31);
        uint32_t absd = (uint32_t)((d32 ^ (int32_t)sign) - (int32_t)sign);
        uint32_t nz   = ct_is_nonzero_u32(absd);

        sm2_affine_t Q = T16[0];
        Q.infinity = 0;

        /* absd 取值 1..16；若 absd==0，后面用 nz 控制不更新 R */
        for (uint32_t v = 1; v <= 16; v++) {
            uint32_t match = ct_eq_u32(absd, v);
            sm2_affine_cmov(&Q, &T16[v - 1], match);
        }

        affine_cond_neg_mont(&Q, sign);

        if (!started) {
            if (absd != 0) {
                sm2_affine_mont_to_jacobian_mont(R, &Q);
                started = 1;
            }
        } else {
            sm2_jacobian_t Radd;
            sm2_add_ja_mont(&Radd, R, &Q);
            sm2_jacobian_cmov(R, &Radd, nz);
        }
    }

    if (!started) {
        sm2_jacobian_set_infinity_mont(R);
    }
}

/* wrapper 1:
 * 输入普通 affine 点 P
 * 输出 Montgomery jacobian
 */
void sm2_scalar_mul_window_ct_schemeB_mont(
    sm2_jacobian_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32])
{
    if (sm2_affine_is_infinity(P) || k_is_zero(k)) {
        sm2_jacobian_set_infinity_mont(R);
        return;
    }

    sm2_affine_t Pm;
    sm2_affine_to_mont(&Pm, P);

    sm2_scalar_mul_window_ct_schemeB_mont_core(R, &Pm, k);
}

/* wrapper 2:
 * 输入普通 affine 点 P
 * 输出 affine（沿用你现有接口语义）
 *
 * 这个函数保留，方便原来的 KEX / 其他上层调用不改接口。
 * 但 benchmark 若想公平比较，不要优先测它，因为它包含最后的 affine 恢复。
 */
void sm2_scalar_mul_window_ct_schemeB_mont_to_affine(
    sm2_affine_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32])
{
    sm2_jacobian_t T;
    sm2_affine_t Pm;

    if (sm2_affine_is_infinity(P) || k_is_zero(k)) {
        R->infinity = 1;
        return;
    }

    sm2_affine_to_mont(&Pm, P);
    sm2_scalar_mul_window_ct_schemeB_mont_core(&T, &Pm, k);

    /* 沿用你当前工程已有语义 */
    sm2_jacobian_to_affine_mont(R, &T);
}

/* wrapper 3:
 * 输入 Montgomery affine 点 P_mont
 * 输出 Montgomery jacobian
 *
 * 这是后续 benchmark 最应该直接测的版本：
 * 不做输入转换，也不做 affine 输出恢复。
 */
void sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(
    sm2_jacobian_t* R,
    const sm2_affine_t* P_mont,
    const uint8_t k[32])
{
    sm2_scalar_mul_window_ct_schemeB_mont_core(R, P_mont, k);
}
void sm2_scalar_mul_window_ct_schemeB_mont_jacobian_to_jacobian_mont(
    sm2_jacobian_t*       R,
    const sm2_jacobian_t* P_mont_jac,
    const uint8_t         k[32])
{
    if (sm2_jacobian_is_infinity(P_mont_jac) || k_is_zero(k)) {
        sm2_jacobian_set_infinity_mont(R);
        return;
    }

    /* Jacobian -> affine（仍在 Montgomery 域内，只是做坐标归一化，不出 Mont 域） */
    sm2_affine_t P_mont_aff;
    sm2_jacobian_mont_to_affine_mont(&P_mont_aff, P_mont_jac);

    /* 走原来的 core 路径 */
    sm2_scalar_mul_window_ct_schemeB_mont_core(R, &P_mont_aff, k);
}
#endif /* USE_FP_MONT */