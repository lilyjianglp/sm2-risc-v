/*fp.h*/
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
 *  所有 fp_* 运算的结果都保证规范化到区间 [0, p-1]
 *  其中 p 为 SM2 曲线使用的素数模数。
 * ============================================================ */

#define FP_LIMBS 4

/* 有限域元素类型 */
typedef struct {
    uint64_t v[FP_LIMBS];
} fp_t;

/* ===================== 常量 ===================== */

/* SM2 素域模数 p */
extern const fp_t FP_P;

/* 常用常量 */
extern const fp_t FP_ZERO;   /* 0 */
extern const fp_t FP_ONE;    /* 1 */

/* ===================== 基本工具函数 ===================== */

/* r = a（拷贝） */
void fp_copy(fp_t* r, const fp_t* a);

/* r = 0 */
void fp_set_zero(fp_t* r);

/* r = 1 */
void fp_set_one(fp_t* r);

/* 判断 a 是否为 0（是返回 1，否则返回 0） */
int  fp_is_zero(const fp_t* a);

/* 判断 a 和 b 是否相等（相等返回 1，否则返回 0） */
int  fp_is_equal(const fp_t* a, const fp_t* b);

/* 整数比较（不取模）
 * 返回值：
 *   -1 : a < b
 *    0 : a == b
 *    1 : a > b
 */
int  fp_cmp(const fp_t* a, const fp_t* b);

/* ===================== 有限域运算（模 p） ===================== */

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
/* 对 512 位中间结果 T 进行模约减，得到 256 位结果 r
 * T 通常是 uint64_t T[8] (512位)
 */
void fp_reduce(fp_t* r, const uint64_t* T);
/* r = a^{-1} mod p
 * 注意：a ≠ 0，由调用者保证
 */
void fp_inv(fp_t* r, const fp_t* a);

/* ===================== 编码 / 解码 ===================== */

/* 将 32 字节大端序字节串转换为有限域元素
 * 若输入值 ≥ p，需要做一次模约减
 */
void fp_from_bytes(fp_t* r, const uint8_t in[32]);

/* 将有限域元素转换为 32 字节大端序字节串 */
void fp_to_bytes(uint8_t out[32], const fp_t* a);

#ifdef __cplusplus
}
#endif

#endif /* SM2_FP_H */
