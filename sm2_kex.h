#pragma once
#include <stdint.h>
#include <stddef.h>
#include "sm2_curve.h"

#ifdef __cplusplus
extern "C" {
#endif

/* SM2 KEX 输出的确认值长度（SM3 输出 32 字节） */
#define SM2_KEX_HASH_SIZE 32

/* 点编码：未压缩 0x04 || X(32) || Y(32) */
#define SM2_POINT_ENC_UNCOMPRESSED_SIZE 65

/* ========== 工具：点编码/解码（协议层常用） ========== */
int sm2_point_encode_uncompressed(uint8_t out[SM2_POINT_ENC_UNCOMPRESSED_SIZE],
                                  const sm2_affine_t *P);

int sm2_point_decode_uncompressed(sm2_affine_t *P,
                                  const uint8_t in[SM2_POINT_ENC_UNCOMPRESSED_SIZE]);

/* 检查点是否在 SM2 曲线上（并且非无穷远点） */
int sm2_point_is_valid(const sm2_affine_t *P);

/* 计算 ZA = Hash(ENTL||IDA||a||b||Gx||Gy||xA||yA) */
int sm2_kex_compute_Z(uint8_t Z[SM2_KEX_HASH_SIZE],
                      const uint8_t *id, size_t id_len,
                      const sm2_affine_t *Ppub);

/* ========== 协议：发起方 A ========== */

/* A 生成临时公钥 RA = [rA]G
 * rA 是 32 字节大端标量（建议随机后 mod n，并确保非 0）
 */
int sm2_kex_initiator_gen_RA(sm2_affine_t *RA,
                             const uint8_t rA[32]);

/* A 在收到 RB 后：
 * - 计算共享点 V
 * - 导出密钥 K（klen 字节）
 * - 生成 S1（发给 B）
 * - 并可用对方发回的 S2 做校验
 *
 * 输入：
 *   dA：A 长期私钥（32B big-endian）
 *   PA：A 长期公钥
 *   idA/idB：双方身份
 *   rA/RA：A 临时标量与点
 *   PB：B 长期公钥
 *   RB：B 临时公钥（收到的）
 * 输出：
 *   K：共享密钥
 *   S1：A->B 的确认值
 *   SA：A 侧计算出的“本端确认值”（等价于标准里的 S_A / S2）
 * 校验：
 *   peer_S2 != NULL 时，会验证；验证失败返回 0
 */
int sm2_kex_initiator_compute_key(uint8_t *K, size_t klen,
                                  uint8_t S1[32],
                                  uint8_t SA[32],
                                  const uint8_t *peer_S2, /* 可为 NULL */
                                  const uint8_t dA[32],
                                  const sm2_affine_t *PA,
                                  const uint8_t *idA, size_t idA_len,
                                  const uint8_t rA[32],
                                  const sm2_affine_t *RA,
                                  const sm2_affine_t *PB,
                                  const uint8_t *idB, size_t idB_len,
                                  const sm2_affine_t *RB);

/* ========== 协议：响应方 B ========== */

/* B 在收到 RA 后：
 * - 生成 RB = [rB]G
 * - 计算共享点 U
 * - 导出密钥 K
 * - 验证 A 发来的 S1
 * - 输出 S2（B->A）
 *
 * 输入：
 *   dB / PB / idB：B 长期私钥、公钥、身份
 *   rB：B 临时标量
 *   RA：收到的 A 临时公钥
 *   PA / idA：A 长期公钥、身份
 *   peer_S1：收到的 S1（必须提供，用于验证）
 * 输出：
 *   RB：B 临时公钥（发回给 A）
 *   K：共享密钥
 *   S2：B->A 的确认值
 *   SB：B 侧计算出的“本端确认值”（等价于标准里的 S_B / S1）
 *
 * 返回：
 *   1 成功（且 S1 验证通过）
 *   0 失败
 */
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

#ifdef __cplusplus
}
#endif

