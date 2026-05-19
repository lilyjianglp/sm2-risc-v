#include "sm2_curve.h"
#include <string.h>

/* ============================================================
 *  基础点函数（普通域，供 sm2_scalar_mont.c / sm2_kex_mont.c 使用）
 *  这些函数只依赖 fp_set_zero/fp_set_one/fp_copy/fp_is_zero，
 *  不调用任何普通域 fp_mul/fp_sqr/fp_inv。
 * ============================================================ */

int sm2_affine_is_infinity(const sm2_affine_t* p) {
    return p->infinity != 0;
}

int sm2_jacobian_is_infinity(const sm2_jacobian_t* p) {
    return fp_is_zero(&p->Z);
}

void sm2_jacobian_set_infinity(sm2_jacobian_t* p) {
    fp_set_zero(&p->X);
    fp_set_one(&p->Y);
    fp_set_zero(&p->Z);
}

/* ============================================================
 *  SM2 曲线常量（Montgomery 域）
 * ============================================================ */
static fp_t SM2_A_MONT;
static fp_t SM2_B_MONT;
static fp_t SM2_GX_MONT;
static fp_t SM2_GY_MONT;

/* 高频小常量，初始化一次后复用，避免每次 double 里重复构造 */
static fp_t SM2_CONST_2_MONT;
static fp_t SM2_CONST_3_MONT;
static fp_t SM2_CONST_8_MONT;

static int g_sm2_curve_mont_inited = 0;


/* 直接复用标准参数字节串 */
static const uint8_t SM2_A_BYTES[32] = {
    0xFF,0xFF,0xFF,0xFE,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0x00,0x00,0x00,0x00,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFC
};

static const uint8_t SM2_B_BYTES[32] = {
    0x28,0xE9,0xFA,0x9E,0x9D,0x9F,0x5E,0x34,
    0x4D,0x5A,0x9E,0x4B,0xCF,0x65,0x09,0xA7,
    0xF3,0x97,0x89,0xF5,0x15,0xAB,0x8F,0x92,
    0xDD,0xBC,0xBD,0x41,0x4D,0x94,0x0E,0x93
};

static const uint8_t SM2_GX_BYTES[32] = {
    0x32,0xC4,0xAE,0x2C,0x1F,0x19,0x81,0x19,
    0x5F,0x99,0x04,0x46,0x6A,0x39,0xC9,0x94,
    0x8F,0xE3,0x0B,0xBF,0xF2,0x66,0x0B,0xE1,
    0x71,0x5A,0x45,0x89,0x33,0x4C,0x74,0xC7
};

static const uint8_t SM2_GY_BYTES[32] = {
    0xBC,0x37,0x36,0xA2,0xF4,0xF6,0x77,0x9C,
    0x59,0xBD,0xCE,0xE3,0x6B,0x69,0x21,0x53,
    0xD0,0xA9,0x87,0x7C,0xC6,0x2A,0x47,0x40,
    0x02,0xDF,0x32,0xE5,0x21,0x39,0xF0,0xA0
};

/* ============================================================
 *  constant-time helpers
 * ============================================================ */

static inline uint32_t sm2_affine_is_infinity_ct_local(const sm2_affine_t* p) {
    return (uint32_t)(p->infinity != 0);
}

static inline uint32_t fp_is_zero_ct_local(const fp_t* a) {
    uint64_t x = a->v[0] | a->v[1] | a->v[2] | a->v[3];
    return (uint32_t)(x == 0);
}

static inline uint32_t sm2_jacobian_is_infinity_ct_local(const sm2_jacobian_t* p) {
    return fp_is_zero_ct_local(&p->Z);
}

static void fp_cmov_local(fp_t* dst, const fp_t* src, uint32_t move) {
    uint64_t mask = 0ULL - (uint64_t)(move & 1U);
    for (int i = 0; i < 4; i++) {
        uint64_t d = dst->v[i];
        uint64_t s = src->v[i];
        dst->v[i] = (d & ~mask) | (s & mask);
    }
}

static void sm2_affine_cmov_local(sm2_affine_t* dst,
                                  const sm2_affine_t* src,
                                  uint32_t move) {
    fp_cmov_local(&dst->x, &src->x, move);
    fp_cmov_local(&dst->y, &src->y, move);
    if (move & 1U) {
        dst->infinity = src->infinity;
    }
}

static void sm2_jacobian_cmov_local(sm2_jacobian_t* dst,
                                    const sm2_jacobian_t* src,
                                    uint32_t move) {
    fp_cmov_local(&dst->X, &src->X, move);
    fp_cmov_local(&dst->Y, &src->Y, move);
    fp_cmov_local(&dst->Z, &src->Z, move);
}

static void sm2_affine_set_infinity_local(sm2_affine_t* p) {
    fp_set_zero(&p->x);
    fp_set_zero(&p->y);
    p->infinity = 1;
}

static void fp_set_small_mont(fp_t* r, uint64_t x) {
    fp_t t;
    fp_set_zero(&t);
    t.v[0] = x;
    fp_to_mont(r, &t);
}

/* ============================================================
 *  初始化
 * ============================================================ */

void sm2_curve_init_once_mont(void) {
    if (g_sm2_curve_mont_inited) return;

    fp_t a, b, gx, gy;

    fp_from_bytes(&a,  SM2_A_BYTES);
    fp_from_bytes(&b,  SM2_B_BYTES);
    fp_from_bytes(&gx, SM2_GX_BYTES);
    fp_from_bytes(&gy, SM2_GY_BYTES);

    fp_to_mont(&SM2_A_MONT,  &a);
    fp_to_mont(&SM2_B_MONT,  &b);
    fp_to_mont(&SM2_GX_MONT, &gx);
    fp_to_mont(&SM2_GY_MONT, &gy);

    /* 这些常量在 doubling/jj 路径里会被高频使用 */
    fp_set_small_mont(&SM2_CONST_2_MONT, 2);
    fp_set_small_mont(&SM2_CONST_3_MONT, 3);
    fp_set_small_mont(&SM2_CONST_8_MONT, 8);

    g_sm2_curve_mont_inited = 1;
}


/* sm2_curve_init_once: alias for mont version */
void sm2_curve_init_once(void) {
    sm2_curve_init_once_mont();
}

/* sm2_get_base_affine: 返回普通域基点（从 Montgomery 域转出）*/
void sm2_get_base_affine(sm2_affine_t* g) {
    sm2_curve_init_once_mont();
    fp_from_mont(&g->x, &SM2_GX_MONT);
    fp_from_mont(&g->y, &SM2_GY_MONT);
    g->infinity = 0;
}

/* ============================================================
 *  基点 / 点表示转换
 * ============================================================ */

void sm2_get_base_affine_mont(sm2_affine_t* g) {
    sm2_curve_init_once_mont();
    fp_copy(&g->x, &SM2_GX_MONT);
    fp_copy(&g->y, &SM2_GY_MONT);
    g->infinity = 0;
}

void sm2_affine_to_mont(sm2_affine_t* r, const sm2_affine_t* p) {
    if (p->infinity) {
        sm2_affine_set_infinity_local(r);
        return;
    }
    fp_to_mont(&r->x, &p->x);
    fp_to_mont(&r->y, &p->y);
    r->infinity = p->infinity;
}

void sm2_affine_from_mont(sm2_affine_t* r, const sm2_affine_t* p) {
    if (p->infinity) {
        sm2_affine_set_infinity_local(r);
        return;
    }
    fp_from_mont(&r->x, &p->x);
    fp_from_mont(&r->y, &p->y);
    r->infinity = p->infinity;
}

void sm2_jacobian_set_infinity_mont(sm2_jacobian_t* p) {
    sm2_curve_init_once_mont();
    fp_set_zero(&p->X);
    fp_copy(&p->Y, &FP_MONT_ONE);
    fp_set_zero(&p->Z);
}

void sm2_affine_to_jacobian_mont(sm2_jacobian_t* r, const sm2_affine_t* q) {
    sm2_curve_init_once_mont();

    fp_to_mont(&r->X, &q->x);
    fp_to_mont(&r->Y, &q->y);

    fp_t zero;
    fp_set_zero(&zero);

    fp_copy(&r->Z, &FP_MONT_ONE);
    fp_cmov_local(&r->Z, &zero, (uint32_t)q->infinity);
}

void sm2_affine_mont_to_jacobian_mont(sm2_jacobian_t* r, const sm2_affine_t* q) {
    sm2_curve_init_once_mont();

    fp_copy(&r->X, &q->x);
    fp_copy(&r->Y, &q->y);

    fp_t zero;
    fp_set_zero(&zero);

    fp_copy(&r->Z, &FP_MONT_ONE);
    fp_cmov_local(&r->Z, &zero, (uint32_t)q->infinity);
}

void sm2_jacobian_mont_to_affine_mont(sm2_affine_t* r, const sm2_jacobian_t* p) {
    sm2_curve_init_once_mont();

    fp_t z, zinv, zinv2, zinv3;
    fp_t x_bar, y_bar;
    sm2_affine_t tmp, inf;

    uint32_t pinf = sm2_jacobian_is_infinity_ct_local(p);

    fp_copy(&z, &p->Z);
    fp_cmov_local(&z, &FP_MONT_ONE, pinf);

    fp_mont_inv(&zinv, &z);
    fp_mont_sqr(&zinv2, &zinv);
    fp_mont_mul(&zinv3, &zinv2, &zinv);

    fp_mont_mul(&x_bar, &p->X, &zinv2);
    fp_mont_mul(&y_bar, &p->Y, &zinv3);

    /* 注意：这里不做 fp_from_mont，保持在 Montgomery 域 */
    fp_copy(&tmp.x, &x_bar);
    fp_copy(&tmp.y, &y_bar);
    tmp.infinity = 0;

    sm2_affine_set_infinity_local(&inf);
    sm2_affine_cmov_local(&tmp, &inf, pinf);

    *r = tmp;
}

void sm2_jacobian_to_affine_mont(sm2_affine_t* r, const sm2_jacobian_t* p) {
    sm2_curve_init_once_mont();

    fp_t z, zinv, zinv2, zinv3;
    fp_t x_bar, y_bar;
    sm2_affine_t tmp, inf;

    uint32_t pinf = sm2_jacobian_is_infinity_ct_local(p);

    fp_copy(&z, &p->Z);
    fp_cmov_local(&z, &FP_MONT_ONE, pinf);

    fp_mont_inv(&zinv, &z);
    fp_mont_sqr(&zinv2, &zinv);
    fp_mont_mul(&zinv3, &zinv2, &zinv);

    fp_mont_mul(&x_bar, &p->X, &zinv2);
    fp_mont_mul(&y_bar, &p->Y, &zinv3);

    fp_from_mont(&tmp.x, &x_bar);
    fp_from_mont(&tmp.y, &y_bar);
    tmp.infinity = 0;

    sm2_affine_set_infinity_local(&inf);
    sm2_affine_cmov_local(&tmp, &inf, pinf);

    *r = tmp;
}

/* ============================================================
 *  曲线运算核（输入输出都在 Montgomery 域中）
 * ============================================================ */

void sm2_double_jm_a_minus3_mont(sm2_jacobian_t* r, const sm2_jacobian_t* p) {
    uint32_t pinf = sm2_jacobian_is_infinity_ct_local(p);

    /* A = X^2 在这套 a=-3 公式里没有被使用，删掉省一次 mont_sqr */
    fp_t B, C, D, E;
    fp_t ZZ, t1, t2;
    fp_t X3, Y3, Z3;

    sm2_curve_init_once_mont();

    fp_mont_sqr(&B, &p->Y);
    fp_mont_sqr(&C, &B);

    /* D = 4 * X * Y^2 */
    fp_mont_mul(&D, &p->X, &B);
    fp_mont_add(&D, &D, &D);
    fp_mont_add(&D, &D, &D);

    fp_mont_sqr(&ZZ, &p->Z);

    /* E = 3 * (X - Z^2) * (X + Z^2) */
    fp_mont_sub(&t1, &p->X, &ZZ);
    fp_mont_add(&t2, &p->X, &ZZ);
    fp_mont_mul(&E, &t1, &t2);
    fp_mont_mul(&E, &E, &SM2_CONST_3_MONT);

    fp_mont_sqr(&X3, &E);
    fp_mont_sub(&X3, &X3, &D);
    fp_mont_sub(&X3, &X3, &D);

    fp_mont_sub(&t1, &D, &X3);
    fp_mont_mul(&Y3, &E, &t1);

    fp_mont_mul(&t2, &C, &SM2_CONST_8_MONT);
    fp_mont_sub(&Y3, &Y3, &t2);

    fp_mont_mul(&Z3, &p->Y, &p->Z);
    fp_mont_mul(&Z3, &Z3, &SM2_CONST_2_MONT);

    r->X = X3;
    r->Y = Y3;
    r->Z = Z3;

    sm2_jacobian_t inf;
    sm2_jacobian_set_infinity_mont(&inf);
    sm2_jacobian_cmov_local(r, &inf, pinf);
}


/*
 * sm2_add_ja_mont 修正版
 *
 * 原版 bug：使用了 add-2007-bl 公式（带 I=4HH、Rv=2*(S2-Y1) 的变体），
 * 与普通域版 sm2_add_ja 的标准公式不一致，导致结果错误。
 *
 * 修正：完全对齐普通版公式，只把 fp_* 替换为 fp_mont_*。
 *
 * 标准 Jacobian 混合加公式（输入 p=Jacobian, q=Affine/Z=1）：
 *   Z1Z1  = Z1^2
 *   U2    = x2 * Z1Z1
 *   Z1Z1Z1= Z1Z1 * Z1
 *   S2    = y2 * Z1Z1Z1
 *   H     = U2 - X1
 *   R     = S2 - Y1          ← 注意：没有 *2
 *   Z3    = Z1 * H
 *   HH    = H^2
 *   HHH   = H^3
 *   X1HH  = X1 * HH
 *   X3    = R^2 - HHH - 2*X1HH
 *   Y3    = R*(X1HH - X3) - Y1*HHH
 */
void sm2_add_ja_mont(sm2_jacobian_t* r,
                     const sm2_jacobian_t* p,
                     const sm2_affine_t* q_mont)
{
    /* 这里先初始化一次，后面直接复用 FP_MONT_ONE 等常量 */
    sm2_curve_init_once_mont();

    uint32_t pinf = sm2_jacobian_is_infinity_ct_local(p);
    uint32_t qinf = sm2_affine_is_infinity_ct_local(q_mont);

    fp_t Z1Z1, Z1Z1Z1;
    fp_t U2, S2;
    fp_t H, Rv;
    fp_t HH, HHH, X1HH, X1HH2;
    fp_t RR;
    fp_t X3, Y3, Z3;
    fp_t t;

    /* --------------------------------------------------------
     * 通用 mixed-add 主路径
     *
     * 输入:
     *   p = Jacobian
     *   q = affine (Montgomery 域)
     *
     * 标准公式:
     *   Z1Z1   = Z1^2
     *   U2     = x2 * Z1^2
     *   Z1Z1Z1 = Z1^3
     *   S2     = y2 * Z1^3
     *   H      = U2 - X1
     *   R      = S2 - Y1
     *   Z3     = Z1 * H
     *   HH     = H^2
     *   HHH    = H^3
     *   X1HH   = X1 * HH
     *   X3     = R^2 - HHH - 2*X1HH
     *   Y3     = R*(X1HH - X3) - Y1*HHH
     * -------------------------------------------------------- */
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

    /* 先写通用加法结果 */
    r->X = X3;
    r->Y = Y3;
    r->Z = Z3;

    /* --------------------------------------------------------
     * 特殊情况掩码
     * -------------------------------------------------------- */
    uint32_t h_zero  = fp_is_zero_ct_local(&H);
    uint32_t r_zero  = fp_is_zero_ct_local(&Rv);
    uint32_t use_dbl = h_zero & r_zero;
    uint32_t use_inf = h_zero & (r_zero ^ 1U);

    /* --------------------------------------------------------
     * 候选结果 1: p == q -> dbl
     *
     * 这是当前 CT 写法里最大的额外成本来源。
     * 这里先保留 constant-time 语义，不引入秘密相关分支。
     * -------------------------------------------------------- */
    sm2_jacobian_t Rdbl;
    sm2_double_jm_a_minus3_mont(&Rdbl, p);

    /* --------------------------------------------------------
     * 候选结果 2: 结果为无穷远
     *
     * 直接就地构造，避免额外 helper 调用和局部临时对象开销。
     * -------------------------------------------------------- */
    sm2_jacobian_t Rinf;
    fp_set_zero(&Rinf.X);
    fp_copy(&Rinf.Y, &FP_MONT_ONE);
    fp_set_zero(&Rinf.Z);

    /* --------------------------------------------------------
     * 候选结果 3: p = infinity 时，结果就是 q
     *
     * 直接内联 affine->jacobian(Mont) 转换，少一次函数调用。
     * -------------------------------------------------------- */
    sm2_jacobian_t Rq;
    fp_copy(&Rq.X, &q_mont->x);
    fp_copy(&Rq.Y, &q_mont->y);
    fp_copy(&Rq.Z, &FP_MONT_ONE);
    {
        fp_t zero;
        fp_set_zero(&zero);
        fp_cmov_local(&Rq.Z, &zero, qinf);
    }

    /* --------------------------------------------------------
     * 按优先级覆盖
     *
     * 默认: 通用 add
     * 然后:
     *   use_dbl -> Rdbl
     *   use_inf -> infinity
     *   pinf    -> q
     *   qinf    -> p
     * -------------------------------------------------------- */
    sm2_jacobian_cmov_local(r, &Rdbl, use_dbl);
    sm2_jacobian_cmov_local(r, &Rinf, use_inf);
    sm2_jacobian_cmov_local(r, &Rq,   pinf);
    sm2_jacobian_cmov_local(r, p,     qinf);
}

void sm2_add_jj_mont(sm2_jacobian_t* r, const sm2_jacobian_t* p, const sm2_jacobian_t* q) {
    uint32_t pinf = sm2_jacobian_is_infinity_ct_local(p);
    uint32_t qinf = sm2_jacobian_is_infinity_ct_local(q);

    fp_t Z1Z1, Z2Z2, U1, U2, S1, S2;
    fp_t H, I, J, Rv, V;
    fp_t X3, Y3, Z3;

    sm2_curve_init_once_mont();

    fp_mont_sqr(&Z1Z1, &p->Z);
    fp_mont_sqr(&Z2Z2, &q->Z);

    fp_mont_mul(&U1, &p->X, &Z2Z2);
    fp_mont_mul(&U2, &q->X, &Z1Z1);

    fp_t Z2_cu, Z1_cu;
    fp_mont_mul(&Z2_cu, &Z2Z2, &q->Z);
    fp_mont_mul(&Z1_cu, &Z1Z1, &p->Z);

    fp_mont_mul(&S1, &p->Y, &Z2_cu);
    fp_mont_mul(&S2, &q->Y, &Z1_cu);

    fp_mont_sub(&H, &U2, &U1);
    fp_mont_add(&I, &H, &H);
    fp_mont_sqr(&I, &I);
    fp_mont_mul(&J, &H, &I);

    fp_mont_sub(&Rv, &S2, &S1);
    fp_mont_add(&Rv, &Rv, &Rv);

    fp_mont_mul(&V, &U1, &I);

    fp_mont_sqr(&X3, &Rv);
    fp_mont_sub(&X3, &X3, &J);
    fp_mont_sub(&X3, &X3, &V);
    fp_mont_sub(&X3, &X3, &V);

    fp_mont_sub(&Y3, &V, &X3);
    fp_mont_mul(&Y3, &Y3, &Rv);

    fp_t S1J2;
    fp_mont_mul(&S1J2, &S1, &J);
    fp_mont_add(&S1J2, &S1J2, &S1J2);
    fp_mont_sub(&Y3, &Y3, &S1J2);

    fp_mont_add(&Z3, &p->Z, &q->Z);
    fp_mont_sqr(&Z3, &Z3);
    fp_mont_sub(&Z3, &Z3, &Z1Z1);
    fp_mont_sub(&Z3, &Z3, &Z2Z2);
    fp_mont_mul(&Z3, &Z3, &H);

    r->X = X3;
    r->Y = Y3;
    r->Z = Z3;

    sm2_jacobian_cmov_local(r, q, pinf);
    sm2_jacobian_cmov_local(r, p, qinf);

    uint32_t h_zero = fp_is_zero_ct_local(&H);
    fp_t SS;
    fp_mont_sub(&SS, &S2, &S1);
    uint32_t s_zero = fp_is_zero_ct_local(&SS);

    uint32_t not_pinf = pinf ^ 1U;
    uint32_t not_qinf = qinf ^ 1U;
    uint32_t not_s_zero = s_zero ^ 1U;

    sm2_jacobian_t dbl, inf;
    sm2_double_jm_a_minus3_mont(&dbl, p);
    sm2_jacobian_set_infinity_mont(&inf);

    sm2_jacobian_cmov_local(r, &dbl, h_zero & s_zero & not_pinf & not_qinf);
    sm2_jacobian_cmov_local(r, &inf, h_zero & not_s_zero & not_pinf & not_qinf);
}


void sm2_batch_normalize_mont(sm2_affine_t* out,
                              const sm2_jacobian_t* in,
                              size_t n,
                              fp_t* prefix_buf) {
    sm2_curve_init_once_mont();

    if (n == 0) return;

    fp_t acc;
    fp_copy(&acc, &FP_MONT_ONE);

    for (size_t i = 0; i < n; i++) {
        prefix_buf[i] = acc;

        uint32_t inf = sm2_jacobian_is_infinity_ct_local(&in[i]);

        fp_t z_or_one;
        fp_copy(&z_or_one, &in[i].Z);
        fp_cmov_local(&z_or_one, &FP_MONT_ONE, inf);

        fp_mont_mul(&acc, &acc, &z_or_one);
    }

    fp_t acc_inv;
    fp_mont_inv(&acc_inv, &acc);

    for (size_t i = n; i-- > 0; ) {
        uint32_t inf = sm2_jacobian_is_infinity_ct_local(&in[i]);

        fp_t z_or_one;
        fp_copy(&z_or_one, &in[i].Z);
        fp_cmov_local(&z_or_one, &FP_MONT_ONE, inf);

        fp_t zinv;
        fp_mont_mul(&zinv, &acc_inv, &prefix_buf[i]);

        fp_mont_mul(&acc_inv, &acc_inv, &z_or_one);

        fp_t zinv2, zinv3, x_bar, y_bar;
        fp_mont_sqr(&zinv2, &zinv);
        fp_mont_mul(&zinv3, &zinv2, &zinv);

        fp_mont_mul(&x_bar, &in[i].X, &zinv2);
        fp_mont_mul(&y_bar, &in[i].Y, &zinv3);

        fp_copy(&out[i].x, &x_bar);
        fp_copy(&out[i].y, &y_bar);
        out[i].infinity = 0;

        sm2_affine_t inf_pt;
        sm2_affine_set_infinity_local(&inf_pt);
        sm2_affine_cmov_local(&out[i], &inf_pt, inf);
    }
}
