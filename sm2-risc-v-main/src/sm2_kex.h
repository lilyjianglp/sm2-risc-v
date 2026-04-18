#ifndef SM2_KEX_H
#define SM2_KEX_H

#include <stddef.h>
#include <stdint.h>
#include "sm2_curve.h"

#define SM2_POINT_ENC_UNCOMPRESSED_SIZE 65

/* ================= 普通版接口 ================= */

int sm2_kex_initiator_gen_RA(sm2_affine_t *RA,
                             const uint8_t rA[32]);

int sm2_kex_initiator_compute_key(uint8_t *K, size_t klen,
                                  uint8_t S1[32],
                                  uint8_t SA[32],
                                  const uint8_t *peer_S2,
                                  const uint8_t dA[32],
                                  const sm2_affine_t *PA,
                                  const uint8_t *idA, size_t idA_len,
                                  const uint8_t rA[32],
                                  const sm2_affine_t *RA,
                                  const sm2_affine_t *PB,
                                  const uint8_t *idB, size_t idB_len,
                                  const sm2_affine_t *RB);

int sm2_kex_responder_compute_key(sm2_affine_t *RB,
                                  uint8_t *K, size_t klen,
                                  uint8_t S2[32],
                                  uint8_t SB[32],
                                  const uint8_t peer_S1[32],
                                  const uint8_t dB[32],
                                  const sm2_affine_t *PB,
                                  const uint8_t *idB, size_t idB_len,
                                  const uint8_t rB[32],
                                  const sm2_affine_t *RA,
                                  const sm2_affine_t *PA,
                                  const uint8_t *idA, size_t idA_len);
/* 接受普通域 Jacobian 输入，输出普通域 Jacobian
 * 用于 KEX 中间步骤，避免不必要的 jacobian->affine 往返 */
 // 在 sm2_scalar.h 中
int k_is_zero(const uint8_t k[32]);
void sm2_scalar_mul_window_ct_schemeB_jacobian_to_jacobian(
    sm2_jacobian_t*       R,
    const sm2_jacobian_t* P_jac,
    const uint8_t         k[32]);

/* ================= Montgomery版接口 ================= */

int sm2_kex_initiator_gen_RA_mont(sm2_affine_t *RA,
                                  const uint8_t rA[32]);

int sm2_kex_initiator_compute_key_mont(uint8_t *K, size_t klen,
                                       uint8_t S1[32],
                                       uint8_t SA[32],
                                       const uint8_t *peer_S2,
                                       const uint8_t dA[32],
                                       const sm2_affine_t *PA,
                                       const uint8_t *idA, size_t idA_len,
                                       const uint8_t rA[32],
                                       const sm2_affine_t *RA,
                                       const sm2_affine_t *PB,
                                       const uint8_t *idB, size_t idB_len,
                                       const sm2_affine_t *RB);

int sm2_kex_responder_compute_key_mont(sm2_affine_t *RB,
                                       uint8_t *K, size_t klen,
                                       uint8_t S2[32],
                                       uint8_t SB[32],
                                       const uint8_t peer_S1[32],
                                       const uint8_t dB[32],
                                       const sm2_affine_t *PB,
                                       const uint8_t *idB, size_t idB_len,
                                       const uint8_t rB[32],
                                       const sm2_affine_t *RA,
                                       const sm2_affine_t *PA,
                                       const uint8_t *idA, size_t idA_len);

#endif