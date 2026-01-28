#pragma once
#include <stdint.h>
#include "sm2_curve.h"

void sm2_scalar_mul_naf(sm2_jacobian_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32]);

void sm2_scalar_mul_ladder_ct(sm2_jacobian_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32]);
void sm2_scalar_mul_window_ct_schemeB(sm2_jacobian_t* R,
    const sm2_affine_t* P,
    const uint8_t k[32]);