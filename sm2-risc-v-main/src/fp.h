#ifndef SM2_FP_H
#define SM2_FP_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 *  SM2 素域 Fp 的基本运算接口
 *
 *  数值表示方式：
 *    使用 4 个 64 位无符号整数表示一个 256 位大整数
 *    采用小端 limb 存储：
 *
 *      a = v[0]
 *        + 2^64  * v[1]
 *        + 2^128 * v[2]
 *        + 2^192 * v[3]
 *
 *  注意：
 *  本头文件同时提供两套接口：
 *
 *  1) 普通域接口（standard domain）
 *     元素表示普通的 a mod p
 *
 *  2) Montgomery 域接口（Montgomery domain）
 *     元素表示 aR mod p，其中 R = 2^256
 *
 *  推荐约定：
 *    - 对外输入/输出、编码/解码：使用普通域接口
 *    - 曲线内部点运算：使用 Montgomery 域接口
 * ============================================================ */

#define FP_LIMBS 4

typedef struct {
    uint64_t v[FP_LIMBS];
} fp_t;

/* ===================== 常量 ===================== */

/* SM2 素域模数 p */
extern const fp_t FP_P;

/* 普通域常量 */
extern const fp_t FP_ZERO;   /* 0 */
extern const fp_t FP_ONE;    /* 1 */

/*
 * Montgomery 域常量：
 * FP_MONT_ONE = R mod p，对应普通域中的 1
 * FP_MONT_R2  = R^2 mod p，用于 to_mont 转换
 */
extern const fp_t FP_MONT_ONE;
extern const fp_t FP_MONT_R2;

/* ===================== 基本工具函数 ===================== */

/* r = a（拷贝） */
void fp_copy(fp_t* r, const fp_t* a);

/* 普通域置 0 / 1 */
void fp_set_zero(fp_t* r);
void fp_set_one(fp_t* r);

/* Montgomery 域置 0 / 1 */
void fp_mont_set_zero(fp_t* r);
void fp_mont_set_one(fp_t* r);

/* 判断是否为 0（对普通域 / Montgomery 域都成立） */
int fp_is_zero(const fp_t* a);

/* 判断是否相等（要求 a,b 在同一表示域内） */
int fp_is_equal(const fp_t* a, const fp_t* b);

/* 整数比较（不取模；要求同一表示域） */
int fp_cmp(const fp_t* a, const fp_t* b);

/* ===================== 普通域运算（模 p） ===================== */

/*
 * 以下接口的输入输出都解释为普通域元素 a mod p
 * 主要用于：
 *   - 编码/解码边界
 *   - 与旧实现/参考实现对拍
 *   - 必要时的外部普通域操作
 */

/* r = (a + b) mod p */
void fp_add(fp_t* r, const fp_t* a, const fp_t* b);

/* r = (a - b) mod p */
void fp_sub(fp_t* r, const fp_t* a, const fp_t* b);

/* r = (-a) mod p */
void fp_neg(fp_t* r, const fp_t* a);

/* r = (a * b) mod p */
void fp_mul(fp_t* r, const fp_t* a, const fp_t* b);

/* r = (a * a) mod p */
void fp_sqr(fp_t* r, const fp_t* a);

/* 对 512 位中间结果 T 做普通模约减，得到 r = T mod p */
void fp_reduce(fp_t* r, const uint64_t* T);

/* r = a^{-1} mod p，要求 a != 0 */
void fp_inv(fp_t* r, const fp_t* a);

/* ===================== Montgomery 域接口 ===================== */

/*
 * Montgomery 域中元素表示为 aR mod p。
 *
 * 曲线内部推荐长期使用这一套接口，避免每次乘法进出域。
 *
 * 注意：
 *   fp_mont_add / sub / neg 与普通域公式形式相同，
 *   因为 (aR ± bR) mod p = (a ± b)R mod p。
 */

/* 边界转换：普通域 -> Montgomery 域 */
void fp_to_mont(fp_t* r, const fp_t* a);

/* 边界转换：Montgomery 域 -> 普通域 */
void fp_from_mont(fp_t* r, const fp_t* a);

/* Montgomery 域加减取负 */
void fp_mont_add(fp_t* r, const fp_t* a, const fp_t* b);
void fp_mont_sub(fp_t* r, const fp_t* a, const fp_t* b);
void fp_mont_neg(fp_t* r, const fp_t* a);

/* Montgomery 域乘法 / 平方
 * 输入输出都保持在 Montgomery 域中：
 *   a_bar = aR mod p
 *   b_bar = bR mod p
 *   fp_mont_mul(r, a_bar, b_bar) => (ab)R mod p
 */
void fp_mont_mul(fp_t* r, const fp_t* a, const fp_t* b);
void fp_mont_sqr(fp_t* r, const fp_t* a);

/* Montgomery 域求逆
 * 若 a_bar = aR mod p，则输出 r_bar = a^{-1}R mod p
 */
void fp_mont_inv(fp_t* r, const fp_t* a);

/*
 * Montgomery REDC 内核：
 * 输入 t 为 9-limb 中间值，输出 r = t * R^{-1} mod p
 * 这是内部内核接口，主要供 benchmark / test 使用。
 */
void fp_mont_reduce_ct(fp_t* r, uint64_t t[9]);

/* ===================== benchmark / test only ===================== */

/* 4x4 -> 8 limb 乘法内核 */
void fp_mul_4x4_u512_ct(uint64_t out[8], const fp_t* a, const fp_t* b);

/* 原始 Montgomery 乘法内核（保留给测试） */
void fp_mont_mul_ct(fp_t* r, const fp_t* a, const fp_t* b);
/* ===================== 编码 / 解码 ===================== */

/*
 * 将 32 字节大端序字节串转换为普通域元素
 * 若输入值 >= p，需要做一次模约减
 */
void fp_from_bytes(fp_t* r, const uint8_t in[32]);

/* 将普通域元素转换为 32 字节大端序字节串 */
void fp_to_bytes(uint8_t out[32], const fp_t* a);
#ifdef __cplusplus
}
#endif

#endif /* SM2_FP_H */