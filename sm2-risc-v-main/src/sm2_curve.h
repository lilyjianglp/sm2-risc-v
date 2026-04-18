#pragma once
#include <stdint.h>
#include <stddef.h>
#include "fp.h"

/* ===== 点结构 ===== */
typedef struct {
    fp_t x, y;
    int infinity;
} sm2_affine_t;

typedef struct {
    fp_t X, Y, Z;   /* Z = 0 表示无穷远点 */
} sm2_jacobian_t;

/* ===== 初始化与基础接口 ===== */
void sm2_curve_init_once(void);
void sm2_get_base_affine(sm2_affine_t* g);

int  sm2_affine_is_infinity(const sm2_affine_t* p);
int  sm2_jacobian_is_infinity(const sm2_jacobian_t* p);

void sm2_affine_to_jacobian(sm2_jacobian_t* r, const sm2_affine_t* p);
void sm2_jacobian_to_affine(sm2_affine_t* r, const sm2_jacobian_t* p);

/* ===== 曲线点运算：普通域 ===== */
void sm2_double_jm_a_minus3(sm2_jacobian_t* r, const sm2_jacobian_t* p);
void sm2_add_ja(sm2_jacobian_t* r, const sm2_jacobian_t* p, const sm2_affine_t* q);
void sm2_add_jj(sm2_jacobian_t* r, const sm2_jacobian_t* p, const sm2_jacobian_t* q);

void sm2_batch_normalize(sm2_affine_t* out,
                         const sm2_jacobian_t* in,
                         size_t n,
                         fp_t* prefix_buf);

void sm2_jacobian_set_infinity(sm2_jacobian_t* p);

/* ===================== Montgomery 曲线层接口 ===================== */

void sm2_curve_init_once_mont(void);

void sm2_get_base_affine_mont(sm2_affine_t* g);

void sm2_affine_to_mont(sm2_affine_t* r, const sm2_affine_t* p);
void sm2_affine_from_mont(sm2_affine_t* r, const sm2_affine_t* p);

void sm2_jacobian_set_infinity_mont(sm2_jacobian_t* p);

void sm2_affine_to_jacobian_mont(sm2_jacobian_t* r, const sm2_affine_t* p);
void sm2_affine_mont_to_jacobian_mont(sm2_jacobian_t* r, const sm2_affine_t* p);
void sm2_jacobian_to_affine_mont(sm2_affine_t* r, const sm2_jacobian_t* p);

void sm2_double_jm_a_minus3_mont(sm2_jacobian_t* r, const sm2_jacobian_t* p);
void sm2_add_ja_mont(sm2_jacobian_t* r, const sm2_jacobian_t* p, const sm2_affine_t* q_mont);
void sm2_add_jj_mont(sm2_jacobian_t* r, const sm2_jacobian_t* p, const sm2_jacobian_t* q);

void sm2_batch_normalize_mont(sm2_affine_t* out,
                              const sm2_jacobian_t* in,
                              size_t n,
                              fp_t* prefix_buf);

void sm2_jacobian_mont_to_affine_mont(sm2_affine_t* r,
                                      const sm2_jacobian_t* p);