#pragma once
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 *  SM2 Fn (mod n) scalar arithmetic
 *  n = order of base point G
 * ============================================================ */

#define FN_LIMBS 4

typedef struct {
    uint64_t v[FN_LIMBS]; /* little-endian limbs: v[0] is least significant */
} fn_t;

/* Constructors / conversions */
void fn_zero(fn_t *r);
int  fn_is_zero(const fn_t *a);

void fn_from_be(fn_t *r, const uint8_t in[32]);  /* 32-byte big-endian -> fn_t (not automatically reduced) */
void fn_to_be(uint8_t out[32], const fn_t *a);   /* fn_t -> 32-byte big-endian */

int  fn_cmp(const fn_t *a, const fn_t *b);       /* -1/0/1 */

/* Reduce into [0, n) */
void fn_normalize(fn_t *a);

/* Arithmetic mod n */
void fn_add(fn_t *r, const fn_t *a, const fn_t *b);      /* r = (a + b) mod n */
void fn_sub(fn_t *r, const fn_t *a, const fn_t *b);      /* r = (a - b) mod n */
void fn_mul(fn_t *r, const fn_t *a, const fn_t *b);      /* r = (a * b) mod n */
void fn_mul_4x4_u512_ct(uint64_t out[8], const fn_t *a, const fn_t *b);
void fn_mont_reduce_ct(fn_t *r, uint64_t t[9]);
void fn_mont_mul_ct(fn_t *r, const fn_t *a, const fn_t *b);
void fn_to_mont(fn_t *r, const fn_t *a);
void fn_from_mont(fn_t *r, const fn_t *a);

/* 新增：Montgomery 连续运算接口 */
void fn_to_mont_pub(fn_t *r, const fn_t *a);
void fn_from_mont_pub(fn_t *r, const fn_t *a);
void fn_mont_mul_pub(fn_t *r, const fn_t *a, const fn_t *b);
void fn_mont_sqr_pub(fn_t *r, const fn_t *a);
void fn_mont_one(fn_t *r);



/* Convenience: compute t = (d + x * r) mod n, where inputs are 32-byte big-endian.
 * Return 1 on success, 0 if result is 0 (or inputs are 0 after normalization).
 */
int fn_compute_t(uint8_t t_be[32],
                 const uint8_t d_be[32],
                 const uint8_t r_be[32],
                 const fn_t *x);

/* x' for SM2 KEX: w = 127, x' = 2^w + (x mod 2^w).
 * Input x is 32-byte big-endian (usually x-coordinate bytes).
 */
void fn_x_dash_w127(fn_t *x_dash, const uint8_t x_be[32]);

/* Expose constant n (order) if needed */
const fn_t* fn_modulus_n(void);

#ifdef __cplusplus
}
#endif
