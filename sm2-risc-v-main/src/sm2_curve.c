/* sm2_curve.c
 *
 * 单文件实现：坐标层 + 曲线代数层
 * 依赖：fp.h / fp.c 提供 fp_* 接口
 */

#include "sm2_curve.h"
#include <stdlib.h>
#include <string.h>
/* ===== constant-time helpers forward declarations ===== */
static inline uint32_t sm2_jacobian_is_infinity_ct(const sm2_jacobian_t* p);



/* ============================================================
 *  1) 曲线参数（SM2 推荐曲线）
 *     y^2 = x^3 + ax + b over Fp
 *     a = -3 (mod p)
 * ============================================================ */

static fp_t SM2_A;
static fp_t SM2_B;
static fp_t SM2_GX;
static fp_t SM2_GY;

static int g_sm2_curve_inited = 0;

/* 标准参数（大端 32 字节）sm2p256v1 */
static const uint8_t SM2_A_BYTES[32] = {
    /* a = FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC */
    /* a = -3 */
    0xFF,0xFF,0xFF,0xFE,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0x00,0x00,0x00,0x00,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFC
};

static const uint8_t SM2_B_BYTES[32] = {
    /* b = 28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93 */
    0x28,0xE9,0xFA,0x9E,0x9D,0x9F,0x5E,0x34,
    0x4D,0x5A,0x9E,0x4B,0xCF,0x65,0x09,0xA7,
    0xF3,0x97,0x89,0xF5,0x15,0xAB,0x8F,0x92,
    0xDD,0xBC,0xBD,0x41,0x4D,0x94,0x0E,0x93
};
/* 标准基点参数 G(G_x  ,G_y​ )*/
static const uint8_t SM2_GX_BYTES[32] = {
    /* Gx = 32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7 */
    0x32,0xC4,0xAE,0x2C,0x1F,0x19,0x81,0x19,
    0x5F,0x99,0x04,0x46,0x6A,0x39,0xC9,0x94,
    0x8F,0xE3,0x0B,0xBF,0xF2,0x66,0x0B,0xE1,
    0x71,0x5A,0x45,0x89,0x33,0x4C,0x74,0xC7
};

static const uint8_t SM2_GY_BYTES[32] = {
    /* Gy = BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0 */
    0xBC,0x37,0x36,0xA2,0xF4,0xF6,0x77,0x9C,
    0x59,0xBD,0xCE,0xE3,0x6B,0x69,0x21,0x53,
    0xD0,0xA9,0x87,0x7C,0xC6,0x2A,0x47,0x40,
    0x02,0xDF,0x32,0xE5,0x21,0x39,0xF0,0xA0
};
/*初始化函数*/
void sm2_curve_init_once(void)
{
    if (g_sm2_curve_inited) return;
    fp_from_bytes(&SM2_A, SM2_A_BYTES);
    fp_from_bytes(&SM2_B, SM2_B_BYTES);
    fp_from_bytes(&SM2_GX, SM2_GX_BYTES);
    fp_from_bytes(&SM2_GY, SM2_GY_BYTES);
    g_sm2_curve_inited = 1;
}
static inline uint32_t sm2_affine_is_infinity_ct(const sm2_affine_t* p)
{
    return (uint32_t)(p->infinity);
}

/* ============================================================
 *  2) 点结构与工具函数（坐标层）
 * ============================================================ */
static void fp_cmov(fp_t* dst, const fp_t* src, uint32_t move)
{
    uint64_t mask = 0ULL - (uint64_t)(move & 1U);
    for (int i = 0; i < 4; i++) {
        uint64_t d = dst->v[i];
        uint64_t s = src->v[i];
        dst->v[i] = (d & ~mask) | (s & mask);
    }
}

static void sm2_affine_cmov(sm2_affine_t* dst,
                            const sm2_affine_t* src,
                            uint32_t move)
{
    fp_cmov(&dst->x, &src->x, move);
    fp_cmov(&dst->y, &src->y, move);
    dst->infinity =
        (dst->infinity & (int)(~(move & 1U))) |
        (src->infinity & (int)(move & 1U));
}

/* 实现 Jacobian 点的常量时间选择 (CT) */
static void sm2_jacobian_cmov(sm2_jacobian_t* dst,
                              const sm2_jacobian_t* src,
                              uint32_t move)
{
    fp_cmov(&dst->X, &src->X, move);
    fp_cmov(&dst->Y, &src->Y, move);
    fp_cmov(&dst->Z, &src->Z, move);
}
/* 实现 Affine 到 Jacobian 的转换 (CT) */
void sm2_affine_to_jacobian(sm2_jacobian_t* r, const sm2_affine_t* q)
{
    fp_copy(&r->X, &q->x);
    fp_copy(&r->Y, &q->y);
    
    /* 如果 q 是无穷远点，Z 设为 0，否则设为 1 */
    fp_t one, zero;
    fp_set_one(&one);
    fp_set_zero(&zero);
    
    fp_copy(&r->Z, &one);
    fp_cmov(&r->Z, &zero, (uint32_t)q->infinity);
}
/* affine infinity */
void sm2_affine_set_infinity(sm2_affine_t* p)
{
    fp_set_zero(&p->x);
    fp_set_zero(&p->y);
    p->infinity = 1;
}
int sm2_affine_is_infinity(const sm2_affine_t* p)
{
    return p->infinity != 0;
}

/* jacobian infinity */
void sm2_jacobian_set_infinity(sm2_jacobian_t* p)
{
    fp_set_zero(&p->X);
    fp_set_one(&p->Y);
    fp_set_zero(&p->Z);//Z=0代表无穷远点

}
int sm2_jacobian_is_infinity(const sm2_jacobian_t* p)
{
    return fp_is_zero(&p->Z);/*判断Z是否为0*/
}

/* 坐标转换：Affine -> Jacobian (Z=1) */
void sm2_jacobian_to_affine(sm2_affine_t* r,
                               const sm2_jacobian_t* p)
{
    fp_t one, z, zinv, zinv2, zinv3;
    sm2_affine_t tmp, inf;

    fp_set_one(&one);

    /* pinf = 1 if Z==0 */
    uint32_t pinf = sm2_jacobian_is_infinity_ct(p);

    /* z = (pinf ? 1 : Z) */
    fp_copy(&z, &p->Z);
    fp_cmov(&z, &one, pinf);

    /* zinv = z^{-1} */
    fp_inv(&zinv, &z);

    /* zinv^2, zinv^3 */
    fp_sqr(&zinv2, &zinv);
    fp_mul(&zinv3, &zinv2, &zinv);

    /* tmp = affine(p) */
    fp_mul(&tmp.x, &p->X, &zinv2);
    fp_mul(&tmp.y, &p->Y, &zinv3);
    tmp.infinity = 0;

    /* inf = infinity */
    sm2_affine_set_infinity(&inf);

    /* if pinf, r = inf else r = tmp */
    sm2_affine_cmov(&tmp, &inf, pinf);

    *r = tmp;
}


/* 取基点G（Affine） */
void sm2_get_base_affine(sm2_affine_t* g)
{
    sm2_curve_init_once();
    fp_copy(&g->x, &SM2_GX);
    fp_copy(&g->y, &SM2_GY);
    g->infinity = 0;
}

/* ============================================================
 *  3) 曲线代数层：点运算核
 *     - 倍加：JM(a=-3) 优化（无逆）
 *     - 点加：JA 混合加法（无逆）
 *     - Jacobian+Jacobian 加法（无逆）
 * ============================================================ */

/* r = 2p, Jacobian, a=-3 优化
 * y2=x3−3x+b
 * 推导公式：
 *  A = X^2
 *  B = Y^2
 *  C = B^2
 *  D = 4*X*B
 *  E = 3*(X - Z^2)*(X + Z^2)   // a=-3 trick: 3*(X^2 - Z^4)
 *  X3 = E^2 - 2D               //结果
 *  Y3 = E*(D - X3) - 8C
 *  Z3 = 2*Y*Z
 */
void sm2_double_jm_a_minus3(sm2_jacobian_t* r, const sm2_jacobian_t* p)
{
    /* 移除 if (...) return; 逻辑，改为计算掩码 */
    uint32_t pinf = sm2_jacobian_is_infinity_ct(p);

    fp_t A, B, C, D, E;
    fp_t ZZ, t1, t2;
    fp_t X3, Y3, Z3;

    /* ===== 以下是标准的倍加计算公式，无论 pinf 是否为真都会执行 ===== */
    fp_sqr(&A, &p->X);         /* A = X^2 */
    fp_sqr(&B, &p->Y);         /* B = Y^2 */
    fp_sqr(&C, &B);            /* C = Y^4 */

    fp_mul(&D, &p->X, &B);     /* D = X*Y^2 */
    fp_add(&D, &D, &D);        /* 2D */
    fp_add(&D, &D, &D);        /* D=4X*Y^2 */

    fp_sqr(&ZZ, &p->Z);        /* ZZ = Z^2 */
    fp_sub(&t1, &p->X, &ZZ);   /* t1 = X - Z^2 */
    fp_add(&t2, &p->X, &ZZ);   /* t2 = X + Z^2 */
    fp_mul(&E, &t1, &t2);      /* E = (X-Z^2)(X+Z^2) */

    /* E = 3E */
    fp_add(&t1, &E, &E);       /* 2E */
    fp_add(&E, &t1, &E);       /* 3E */

    fp_sqr(&X3, &E);           /* X3 = E^2 */
    fp_add(&t1, &D, &D);       /* t1 = 2D */
    fp_sub(&X3, &X3, &t1);     /* X3 = E^2 - 2D */

    fp_sub(&t1, &D, &X3);      /* t1 = D - X3 */
    fp_mul(&Y3, &E, &t1);      /* Y3 = E*(D - X3) */

    /* subtract 8C */
    fp_add(&t1, &C, &C);       /* 2C */
    fp_add(&t1, &t1, &t1);     /* 4C */
    fp_add(&t1, &t1, &t1);     /* 8C */
    fp_sub(&Y3, &Y3, &t1);     /* Y3 -= 8C */

    fp_mul(&Z3, &p->Y, &p->Z); /* Z3 = Y*Z */
    fp_add(&Z3, &Z3, &Z3);     /* Z3 = 2YZ */

    /* ===== 结果合并：使用常量时间选择 ===== */
    
    /* 1. 先将计算结果写入 r */
    fp_copy(&r->X, &X3);
    fp_copy(&r->Y, &Y3);
    fp_copy(&r->Z, &Z3);

    /* 2. 如果输入本身是无穷远点，则结果仍应为无穷远点 */
    /* 我们需要一个临时的无穷远点结构 */
    sm2_jacobian_t Tinf;
    sm2_jacobian_set_infinity(&Tinf); 

    /* 使用 cmov：如果 pinf 为 1，则将 r 覆盖为 Tinf */
    sm2_jacobian_cmov(r, &Tinf, pinf);
}

static inline uint32_t ct_is_zero_u64(uint64_t x)
{
    return (uint32_t)(((x | (0ULL - x)) >> 63) ^ 1ULL);
}

static inline uint32_t fp_is_zero_ct(const fp_t* a)
{
    uint64_t acc = 0;
    for (int i = 0; i < 4; i++) acc |= a->v[i];
    return ct_is_zero_u64(acc); /* 1 if zero */
}

static inline uint32_t sm2_jacobian_is_infinity_ct(const sm2_jacobian_t* p)
{
    return fp_is_zero_ct(&p->Z); /* 1 if infinity */
}

/* r = p + q, where p is Jacobian, q is Affine (Z=1)
 *
 * 标准混合加JA：
 *  Z1Z1 = Z1^2
 *  U2 = x2*Z1Z1
 *  Z1Z1Z1 = Z1Z1*Z1
 *  S2 = y2*Z1Z1Z1
 *  H = U2 - X1
 *  R = S2 - Y1
 *  Z3 = Z1*H
 *  HH = H^2
 *  HHH = H^3
 *  X1HH = X1*HH
 *  X3 = R^2 - HHH - 2*X1HH
 *  Y3 = R*(X1HH - X3) - Y1*HHH
 *
 
 */

void sm2_add_ja(sm2_jacobian_t* r,
                const sm2_jacobian_t* p,
                const sm2_affine_t* q)
{
    /* ===== 通用混合加：计算 Radd = p + q（假设 H!=0 且 p 非∞） ===== */
    fp_t Z1Z1, Z1Z1Z1;
    fp_t U2, S2;
    fp_t H, Rv;
    fp_t HH, HHH, X1HH, X1HH2;
    fp_t RR;
    fp_t X3a, Y3a, Z3a;
    fp_t t;

    fp_sqr(&Z1Z1, &p->Z);              /* Z^2 */
    fp_mul(&U2, &q->x, &Z1Z1);         /* x2*Z^2 */

    fp_mul(&Z1Z1Z1, &Z1Z1, &p->Z);     /* Z^3 */
    fp_mul(&S2, &q->y, &Z1Z1Z1);       /* y2*Z^3 */

    fp_sub(&H, &U2, &p->X);            /* H = U2 - X1 */
    fp_sub(&Rv, &S2, &p->Y);           /* R = S2 - Y1 */

    fp_mul(&Z3a, &p->Z, &H);           /* Z3 = Z*H */

    fp_sqr(&HH, &H);                   /* H^2 */
    fp_mul(&HHH, &HH, &H);             /* H^3 */
    fp_mul(&X1HH, &p->X, &HH);         /* X1*H^2 */
    fp_add(&X1HH2, &X1HH, &X1HH);      /* 2*X1HH */

    fp_sqr(&RR, &Rv);                  /* R^2 */

    fp_sub(&t, &RR, &HHH);             /* R^2 - H^3 */
    fp_sub(&X3a, &t, &X1HH2);          /* X3 */

    fp_sub(&t, &X1HH, &X3a);           /* X1HH - X3 */
    fp_mul(&t, &Rv, &t);               /* R*(X1HH - X3) */
    fp_mul(&Y3a, &p->Y, &HHH);         /* Y1*H^3 */
    fp_sub(&Y3a, &t, &Y3a);            /* Y3 */

    /* ===== 候选：Rdbl（H==0 && R==0） ===== */
    sm2_jacobian_t Rdbl;
    sm2_double_jm_a_minus3(&Rdbl, p);

    /* ===== 候选：Rinf（H==0 && R!=0） ===== */
    sm2_jacobian_t Rinf;
    sm2_jacobian_set_infinity(&Rinf);

    /* ===== 候选：Rq（p==∞ 时结果就是 q） ===== */
    sm2_jacobian_t Rq;
    sm2_affine_to_jacobian(&Rq, q);

    /* ===== 条件（全部 constant-time） ===== */
    uint32_t h0   = fp_is_zero_ct(&H);                 /* H==0 */
    uint32_t r0   = fp_is_zero_ct(&Rv);                /* R==0 */
    uint32_t pinf = sm2_jacobian_is_infinity_ct(p);    /* p==∞ */

    uint32_t use_dbl = (uint32_t)(h0 & r0);            /* H==0 && R==0 */
    uint32_t use_inf = (uint32_t)(h0 & (r0 ^ 1U));     /* H==0 && R!=0 */
    /* 默认是 add（H!=0），无需 use_add 掩码 */

    /* ===== 先写入 add 结果 ===== */
    fp_copy(&r->X, &X3a);
    fp_copy(&r->Y, &Y3a);
    fp_copy(&r->Z, &Z3a);

    /* ===== 覆盖：dbl / inf / pinf（顺序很关键） ===== */
    sm2_jacobian_cmov(r, &Rdbl, use_dbl);
    sm2_jacobian_cmov(r, &Rinf, use_inf);
    sm2_jacobian_cmov(r, &Rq,   pinf);
}

/* r = p + q, Jacobian + Jacobian (常量时间版本)
 *
 * 逻辑：
 * 1. 计算通用加法结果 Radd
 * 2. 计算倍加结果 Rdbl (当 p == q)
 * 3. 处理无穷远点情况 (p=inf, q=inf, 或 p=-q)
 * 4. 使用 cmov 选出最终结果
 */
void sm2_add_jj(sm2_jacobian_t* r, const sm2_jacobian_t* p, const sm2_jacobian_t* q)
{
    fp_t Z1Z1, Z2Z2, Z1Z1Z1, Z2Z2Z2;
    fp_t U1, U2, S1, S2;
    fp_t H, Rv, HH, HHH, U1HH;
    fp_t X3, Y3, Z3;
    fp_t t;

    /* ===== 1. 计算中间变量 (无论如何都会执行) ===== */
    fp_sqr(&Z1Z1, &p->Z);          /* Z1^2 mod p */
    fp_sqr(&Z2Z2, &q->Z);          /* Z2^2 mod p */

    fp_mul(&U1, &p->X, &Z2Z2);     /* U1 = X1*Z2^2 */
    fp_mul(&U2, &q->X, &Z1Z1);     /* U2 = X2*Z1^2 */

    fp_mul(&Z1Z1Z1, &Z1Z1, &p->Z); /* Z1^3 */
    fp_mul(&Z2Z2Z2, &Z2Z2, &q->Z); /* Z2^3 */

    fp_mul(&S1, &p->Y, &Z2Z2Z2);   /* S1 = Y1*Z2^3 */
    fp_mul(&S2, &q->Y, &Z1Z1Z1);   /* S2 = Y2*Z1^3 */

    fp_sub(&H, &U2, &U1);          /* H = U2 - U1 */
    fp_sub(&Rv, &S2, &S1);         /* R = S2 - S1 */

    /* ===== 2. 计算通用加法结果 (Radd) ===== */
    fp_sqr(&HH, &H);               /* HH = H^2 */
    fp_mul(&HHH, &HH, &H);         /* HHH = H^3 */
    fp_mul(&U1HH, &U1, &HH);       /* U1HH = U1*H^2 */

    fp_mul(&Z3, &p->Z, &q->Z);
    fp_mul(&Z3, &Z3, &H);          /* Z3 = H*Z1*Z2 */

    fp_sqr(&X3, &Rv);              /* R^2 */
    fp_sub(&X3, &X3, &HHH);
    fp_add(&t, &U1HH, &U1HH);      /* 2*U1HH */
    fp_sub(&X3, &X3, &t);          /* X3 = R^2 - H^3 - 2*U1HH */

    fp_sub(&t, &U1HH, &X3);
    fp_mul(&t, &t, &Rv);
    fp_mul(&Y3, &S1, &HHH);
    fp_sub(&Y3, &t, &Y3);          /* Y3 = R*(U1HH - X3) - S1*H^3 */

    /* ===== 3. 计算备选结果 (候选分支) ===== */
    sm2_jacobian_t Rdbl;
    sm2_double_jm_a_minus3(&Rdbl, p); /* p==q 倍加结果 */

    sm2_jacobian_t Rinf;
    sm2_jacobian_set_infinity(&Rinf); /* p==-q 无穷远结果 */

    /* ===== 4. 判定条件 (Constant-Time) ===== */
    uint32_t pinf = sm2_jacobian_is_infinity_ct(p);
    uint32_t qinf = sm2_jacobian_is_infinity_ct(q);
    uint32_t h0   = fp_is_zero_ct(&H);
    uint32_t r0   = fp_is_zero_ct(&Rv);

    /* ===== 5. 结果筛选 (逻辑优先级由低到高覆盖) ===== */
    // 默认结果设为通用加法结果
    fp_copy(&r->X, &X3);
    fp_copy(&r->Y, &Y3);
    fp_copy(&r->Z, &Z3);

    // 如果 H==0 且 R==0，则是倍加
    sm2_jacobian_cmov(r, &Rdbl, (uint32_t)(h0 & r0));

    // 如果 H==0 且 R!=0，则是互为相反点，结果为无穷远
    sm2_jacobian_cmov(r, &Rinf, (uint32_t)(h0 & (r0 ^ 1U)));

    // 如果 p 是无穷远，结果是 q
    sm2_jacobian_cmov(r, q, pinf);

    // 如果 q 是无穷远，结果是 p (优先级最高)
    sm2_jacobian_cmov(r, p, qinf);
}

/* ============================================================
 *  4) 批量归一化（只在边界用 fp_inv）
 * ============================================================ */

/* out[i] = affine(in[i])，n 个点仅一次逆元
 * 注意：如果 in[i] 是 infinity，out[i] 也置 infinity，并且该点不参与乘积。
 */
/* prefix_buf 长度至少 n+1
 * prefix_buf[0] = 1
 * prefix_buf[i+1] = Π_{j=0..i} (Zj or 1 if infinity)
 * out[i] = affine(in[i])  (CT, no malloc, sentinel avoids idx==0 branch)
 */
void sm2_batch_normalize(sm2_affine_t* out,
                             const sm2_jacobian_t* in,
                             size_t n,
                             fp_t* prefix_buf)
{
    if (n == 0) return;

    fp_t one, acc, t_z;
    fp_set_one(&one);

    /* 哨兵 */
    fp_copy(&prefix_buf[0], &one);

    fp_copy(&acc, &one);

    /* 1) forward: prefix */
    for (size_t i = 0; i < n; i++) {
        uint32_t pinf = sm2_jacobian_is_infinity_ct(&in[i]);

        /* t_z = (pinf ? 1 : Zi) */
        fp_copy(&t_z, &in[i].Z);
        fp_cmov(&t_z, &one, pinf);

        fp_mul(&acc, &acc, &t_z);
        fp_copy(&prefix_buf[i + 1], &acc);
    }

    /* 2) invert total product */
    fp_t inv_acc;
    fp_inv(&inv_acc, &acc); /* 你注释得对：fp_inv 必须是 CT 实现 */

    /* 3) backward */
    for (size_t idx = n; idx-- > 0;) {
        uint32_t pinf = sm2_jacobian_is_infinity_ct(&in[idx]);

        fp_t zinv, zinv2, zinv3;

        /* prev = prefix_buf[idx]  (sentinel: idx=0 => 1) */
        fp_mul(&zinv, &inv_acc, &prefix_buf[idx]);

        /* inv_acc *= (pinf ? 1 : Zi) */
        fp_copy(&t_z, &in[idx].Z);
        fp_cmov(&t_z, &one, pinf);
        fp_mul(&inv_acc, &inv_acc, &t_z);

        /* affine coords */
        fp_sqr(&zinv2, &zinv);
        fp_mul(&zinv3, &zinv2, &zinv);

        sm2_affine_t tmp, inf;
        fp_mul(&tmp.x, &in[idx].X, &zinv2);
        fp_mul(&tmp.y, &in[idx].Y, &zinv3);
        tmp.infinity = 0;

        /* 规范化 infinity 输出（CT 掩蔽） */
        sm2_affine_set_infinity(&inf);
        sm2_affine_cmov(&tmp, &inf, pinf);

        out[idx] = tmp;
    }
}
