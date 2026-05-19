#pragma once
#include <stdint.h>
#include "sm2_curve.h"

/* ===== 普通域标量乘 ===== */
void sm2_scalar_mul_window_ct_schemeB(sm2_jacobian_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32]);

/* ===== Montgomery域标量乘 ===== */

/* 输入普通域 affine P，输出 Montgomery Jacobian */
void sm2_scalar_mul_window_ct_schemeB_mont(sm2_jacobian_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32]);

/* 输入普通域 affine P，输出普通域 affine */
void sm2_scalar_mul_window_ct_schemeB_mont_to_affine(sm2_affine_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32]);

/* 输入 Montgomery affine P，输出 Montgomery Jacobian（benchmark 直接用） */
void sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(
    sm2_jacobian_t* R,
    const sm2_affine_t* P_mont,
    const uint8_t k[32]);

/* 输入 Montgomery Jacobian P，输出 Montgomery Jacobian（KEX 中间步骤用）*/
void sm2_scalar_mul_window_ct_schemeB_mont_jacobian_to_jacobian_mont(
    sm2_jacobian_t*       R,
    const sm2_jacobian_t* P_mont_jac,
    const uint8_t         k[32]);