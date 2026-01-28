#pragma once
#include <stdint.h>
#include <stddef.h>
#include "fp.h"

/* ===== ��ṹ ===== */
typedef struct {
    fp_t x, y;
    int infinity;
} sm2_affine_t;

typedef struct {
    fp_t X, Y, Z;   /* Z=0 => infinity */
} sm2_jacobian_t;

/* ===== ���� & ���� ===== */
void sm2_curve_init_once(void);
void sm2_get_base_affine(sm2_affine_t* g);

int  sm2_affine_is_infinity(const sm2_affine_t* p);
int  sm2_jacobian_is_infinity(const sm2_jacobian_t* p);

void sm2_affine_to_jacobian(sm2_jacobian_t* r, const sm2_affine_t* p);
void sm2_jacobian_to_affine(sm2_affine_t* r, const sm2_jacobian_t* p);

/* ===== �ڶ��㣺������� ===== */
void sm2_double_jm_a_minus3(sm2_jacobian_t* r, const sm2_jacobian_t* p);
void sm2_add_ja(sm2_jacobian_t* r, const sm2_jacobian_t* p, const sm2_affine_t* q);
void sm2_add_jj(sm2_jacobian_t* r, const sm2_jacobian_t* p, const sm2_jacobian_t* q);

void sm2_batch_normalize(sm2_affine_t* out,
                         const sm2_jacobian_t* in,
                         size_t n,
                         fp_t* prefix_buf);

void sm2_jacobian_set_infinity(sm2_jacobian_t* p);