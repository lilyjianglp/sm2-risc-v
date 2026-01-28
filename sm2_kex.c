#include "sm2_kex.h"
#include "fp.h"
#include "sm2_scalar.h"
#include <string.h>
#include "fn.h"

/* ============================================================
 *  SM2 KEX: 内置 SM3 + KDF（不依赖外部库）
 * ============================================================ */

typedef struct {
    uint32_t h[8];
    uint64_t nbits;
    uint8_t  buf[64];
    size_t   buf_len;
} sm3_ctx;

static uint32_t rol32(uint32_t x, uint32_t n) {
    return (x << n) | (x >> (32 - n));
}

static uint32_t load_be32(const uint8_t *p) {
    return ((uint32_t)p[0] << 24) | ((uint32_t)p[1] << 16) | ((uint32_t)p[2] << 8) | (uint32_t)p[3];
}

static void store_be32(uint8_t *p, uint32_t x) {
    p[0] = (uint8_t)(x >> 24);
    p[1] = (uint8_t)(x >> 16);
    p[2] = (uint8_t)(x >> 8);
    p[3] = (uint8_t)x;
}

static uint32_t sm3_p0(uint32_t x) { return x ^ rol32(x, 9) ^ rol32(x, 17); }
static uint32_t sm3_p1(uint32_t x) { return x ^ rol32(x, 15) ^ rol32(x, 23); }

static uint32_t sm3_ff(uint32_t x, uint32_t y, uint32_t z, int j) {
    return (j < 16) ? (x ^ y ^ z) : ((x & y) | (x & z) | (y & z));
}
static uint32_t sm3_gg(uint32_t x, uint32_t y, uint32_t z, int j) {
    return (j < 16) ? (x ^ y ^ z) : ((x & y) | ((~x) & z));
}

static void sm3_compress(uint32_t H[8], const uint8_t block[64]) {
    uint32_t W[68];
    uint32_t W1[64];

    for (int i = 0; i < 16; i++) {
        W[i] = load_be32(block + 4 * i);
    }
    for (int j = 16; j < 68; j++) {
        uint32_t x = W[j - 16] ^ W[j - 9] ^ rol32(W[j - 3], 15);
        W[j] = sm3_p1(x) ^ rol32(W[j - 13], 7) ^ W[j - 6];
    }
    for (int j = 0; j < 64; j++) {
        W1[j] = W[j] ^ W[j + 4];
    }

    uint32_t A = H[0], B = H[1], C = H[2], D = H[3];
    uint32_t E = H[4], F = H[5], G = H[6], HH = H[7];

    for (int j = 0; j < 64; j++) {
        uint32_t Tj = (j < 16) ? 0x79CC4519u : 0x7A879D8Au;
        uint32_t SS1 = rol32((rol32(A, 12) + E + rol32(Tj, j)) & 0xFFFFFFFFu, 7);
        uint32_t SS2 = SS1 ^ rol32(A, 12);
        uint32_t TT1 = (sm3_ff(A, B, C, j) + D + SS2 + W1[j]) & 0xFFFFFFFFu;
        uint32_t TT2 = (sm3_gg(E, F, G, j) + HH + SS1 + W[j]) & 0xFFFFFFFFu;

        D = C;
        C = rol32(B, 9);
        B = A;
        A = TT1;

        HH = G;
        G = rol32(F, 19);
        F = E;
        E = sm3_p0(TT2);
    }

    H[0] ^= A; H[1] ^= B; H[2] ^= C; H[3] ^= D;
    H[4] ^= E; H[5] ^= F; H[6] ^= G; H[7] ^= HH;
}

static void sm3_init(sm3_ctx *ctx) {
    ctx->h[0] = 0x7380166F; ctx->h[1] = 0x4914B2B9;
    ctx->h[2] = 0x172442D7; ctx->h[3] = 0xDA8A0600;
    ctx->h[4] = 0xA96F30BC; ctx->h[5] = 0x163138AA;
    ctx->h[6] = 0xE38DEE4D; ctx->h[7] = 0xB0FB0E4E;
    ctx->nbits = 0;
    ctx->buf_len = 0;
}

static void sm3_update(sm3_ctx *ctx, const uint8_t *in, size_t inlen) {
    ctx->nbits += (uint64_t)inlen * 8;

    while (inlen > 0) {
        size_t n = 64 - ctx->buf_len;
        if (n > inlen) n = inlen;
        memcpy(ctx->buf + ctx->buf_len, in, n);
        ctx->buf_len += n;
        in += n;
        inlen -= n;

        if (ctx->buf_len == 64) {
            sm3_compress(ctx->h, ctx->buf);
            ctx->buf_len = 0;
        }
    }
}

static void sm3_final(sm3_ctx *ctx, uint8_t out[32]) {
    uint8_t pad[72];
    size_t padlen;

    pad[0] = 0x80;
    size_t zlen = (ctx->buf_len < 56) ? (56 - ctx->buf_len - 1) : (64 + 56 - ctx->buf_len - 1);
    memset(pad + 1, 0, zlen);

    uint64_t nbits = ctx->nbits;
    store_be32(pad + 1 + zlen + 0, (uint32_t)(nbits >> 32));
    store_be32(pad + 1 + zlen + 4, (uint32_t)(nbits & 0xFFFFFFFFu));

    padlen = 1 + zlen + 8;
    sm3_update(ctx, pad, padlen);

    for (int i = 0; i < 8; i++) {
        store_be32(out + 4 * i, ctx->h[i]);
    }
}

static void sm3_hash(uint8_t out[32], const uint8_t *in, size_t inlen) {
    sm3_ctx c;
    sm3_init(&c);
    sm3_update(&c, in, inlen);
    sm3_final(&c, out);
}

/* KDF: SM3-based, output klen bytes */
static void sm3_kdf(uint8_t *out, size_t klen, const uint8_t *Z, size_t Zlen) {
    uint32_t ct = 1;
    uint8_t buf[32];
    uint8_t tmp[4];
    size_t off = 0;

    while (off < klen) {
        tmp[0] = (uint8_t)(ct >> 24);
        tmp[1] = (uint8_t)(ct >> 16);
        tmp[2] = (uint8_t)(ct >> 8);
        tmp[3] = (uint8_t)(ct);

        sm3_ctx c;
        sm3_init(&c);
        sm3_update(&c, Z, Zlen);
        sm3_update(&c, tmp, 4);
        sm3_final(&c, buf);

        size_t n = (klen - off < 32) ? (klen - off) : 32;
        memcpy(out + off, buf, n);
        off += n;
        ct++;
    }
}

/* ============================================================
 *  SM2 KEX: 标量域 n 运算
 *
 *  说明：
 *  - 之前 sm2_kex.c 内部自带了一套 256-bit 标量域实现（bn256 + Barrett reduction）。
 *  - 现在这些已经迁移到 fn.c/fn.h（fn_add/fn_sub/fn_mul/fn_normalize、fn_x_dash_w127、fn_compute_t 等）。
 *  - 因此这里删除旧的 bn256 代码，统一复用 fn.c，避免重复实现与不一致风险。
 * ============================================================ */

/* ============================================================
 *  SM2 曲线参数（用于 Z 与点合法性检查）
 * ============================================================ */

/* a,b,Gx,Gy 按标准固定（大端 32B） */
static const uint8_t SM2_A_BYTES[32] = {
    0xFF,0xFF,0xFF,0xFE,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0x00,0x00,0x00,0x00,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFC
};
static const uint8_t SM2_B_BYTES[32] = {
    0x28,0xE9,0xFA,0x9E,0x9D,0x9F,0x5E,0x34,
    0x4D,0x5A,0x9E,0x4B,0xCF,0x65,0x09,0xA7,
    0xF3,0x97,0x89,0xF5,0x15,0xAB,0x8F,0x92,
    0xDD,0xBC,0xBD,0x41,0x4D,0x94,0x0E,0x93
};
static const uint8_t SM2_GX_BYTES[32] = {
    0x32,0xC4,0xAE,0x2C,0x1F,0x19,0x81,0x19,
    0x5F,0x99,0x04,0x46,0x6A,0x39,0xC9,0x94,
    0x8F,0xE3,0x0B,0xBF,0xF2,0x66,0x0B,0xE1,
    0x71,0x5A,0x45,0x89,0x33,0x4C,0x74,0xC7
};
static const uint8_t SM2_GY_BYTES[32] = {
    0xBC,0x37,0x36,0xA2,0xF4,0xF6,0x77,0x9C,
    0x59,0xBD,0xCE,0xE3,0x6B,0x69,0x21,0x53,
    0xD0,0xA9,0x87,0x7C,0xC6,0x2A,0x47,0x40,
    0x02,0xDF,0x32,0xE5,0x21,0x39,0xF0,0xA0
};

static void sm2_curve_params(fp_t *a, fp_t *b) {
    fp_from_bytes(a, SM2_A_BYTES);
    fp_from_bytes(b, SM2_B_BYTES);
}

/* ============================================================
 *  点编码/解码 + 点合法性检查
 * ============================================================ */

int sm2_point_encode_uncompressed(uint8_t out[SM2_POINT_ENC_UNCOMPRESSED_SIZE],
                                  const sm2_affine_t *P)
{
    if (!P || sm2_affine_is_infinity(P)) return 0;
    out[0] = 0x04;
    fp_to_bytes(out + 1, &P->x);
    fp_to_bytes(out + 33, &P->y);
    return 1;
}

int sm2_point_decode_uncompressed(sm2_affine_t *P,
                                  const uint8_t in[SM2_POINT_ENC_UNCOMPRESSED_SIZE])
{
    if (!P || !in) return 0;
    if (in[0] != 0x04) return 0;
    fp_from_bytes(&P->x, in + 1);
    fp_from_bytes(&P->y, in + 33);
    P->infinity = 0;
    return sm2_point_is_valid(P);
}

int sm2_point_is_valid(const sm2_affine_t *P)
{
    if (!P) return 0;
    if (sm2_affine_is_infinity(P)) return 0;

    fp_t a, b;
    sm2_curve_params(&a, &b);

    /* y^2 ?= x^3 + ax + b */
    fp_t y2, x2, x3, ax, rhs;
    fp_sqr(&y2, &P->y);

    fp_sqr(&x2, &P->x);
    fp_mul(&x3, &x2, &P->x);

    fp_mul(&ax, &a, &P->x);
    fp_add(&rhs, &x3, &ax);
    fp_add(&rhs, &rhs, &b);

    return fp_is_equal(&y2, &rhs);
}

/* ============================================================
 *  ZA/ZB 计算
 * ============================================================ */

static void sm2_u16be(uint8_t out[2], uint16_t x) {
    out[0] = (uint8_t)(x >> 8);
    out[1] = (uint8_t)(x);
}

int sm2_kex_compute_Z(uint8_t Z[32],
                      const uint8_t *id, size_t id_len,
                      const sm2_affine_t *Ppub)
{
    if (!Z || !id || !Ppub) return 0;
    if (!sm2_point_is_valid(Ppub)) return 0;

    /* ENTL = id_len * 8 (bits) */
    uint16_t entl = (uint16_t)(id_len * 8);
    uint8_t entl_be[2];
    sm2_u16be(entl_be, entl);

    /* a,b,Gx,Gy */
    fp_t a, b;
    sm2_curve_params(&a, &b);
    uint8_t a_bytes[32], b_bytes[32];
    fp_to_bytes(a_bytes, &a);
    fp_to_bytes(b_bytes, &b);

    sm3_ctx c;
    sm3_init(&c);
    sm3_update(&c, entl_be, 2);
    sm3_update(&c, id, id_len);
    sm3_update(&c, a_bytes, 32);
    sm3_update(&c, b_bytes, 32);
    sm3_update(&c, SM2_GX_BYTES, 32);
    sm3_update(&c, SM2_GY_BYTES, 32);

    uint8_t x_bytes[32], y_bytes[32];
    fp_to_bytes(x_bytes, &Ppub->x);
    fp_to_bytes(y_bytes, &Ppub->y);
    sm3_update(&c, x_bytes, 32);
    sm3_update(&c, y_bytes, 32);

    sm3_final(&c, Z);
    return 1;
}

/* ============================================================
 *  KEX 核心：x' 计算（w=127）
 *  x' = 2^w + (x & (2^w-1))
 * ============================================================ */

static void sm2_kex_x_dash(fn_t *x_dash, const fp_t *x_fp)
{
    uint8_t xb[32];
    fp_to_bytes(xb, x_fp);

    /* x' = 2^127 + (x mod 2^127) */
    fn_x_dash_w127(x_dash, xb);
}

/* ============================================================
 *  计算共享点：
 *    temp = [x_other'] * R_other + P_other
 *    V/U = [t_self] * temp
 *  (SM2 cofactor h=1)
 * ============================================================ */

static int sm2_compute_shared_point(sm2_affine_t *out_xy,
                                   const uint8_t t_self_be[32],
                                   const fn_t *x_other_dash,
                                   const sm2_affine_t *R_other,
                                   const sm2_affine_t *P_other)
{
    if (!out_xy || !t_self_be || !x_other_dash || !R_other || !P_other) return 0;
    if (!sm2_point_is_valid(R_other)) return 0;
    if (!sm2_point_is_valid(P_other)) return 0;

    uint8_t xdash_be[32];
    fn_to_be(xdash_be, x_other_dash);

    /* [x_other'] * R_other */
    sm2_jacobian_t J1;
    sm2_scalar_mul_window_ct_schemeB(&J1, R_other, xdash_be);

    /* + P_other */
    sm2_add_ja(&J1, &J1, P_other);

    /* [t_self] * temp */
    sm2_affine_t temp_aff;
    sm2_jacobian_to_affine(&temp_aff, &J1);
    if (sm2_affine_is_infinity(&temp_aff)) return 0;

    sm2_jacobian_t J2;
    sm2_scalar_mul_window_ct_schemeB(&J2, &temp_aff, t_self_be);

    sm2_jacobian_to_affine(out_xy, &J2);
    if (sm2_affine_is_infinity(out_xy)) return 0;
    return 1;
}

/* ============================================================
 *  计算 t = (d + x_dash * r) mod n
 * ============================================================ */
static int sm2_compute_t(uint8_t t_be[32],
                         const uint8_t d_be[32],
                         const uint8_t r_be[32],
                         const fn_t *x_dash)
{
    return fn_compute_t(t_be, d_be, r_be, x_dash);
}

/* ============================================================
 *  inner = SM3(xV||ZA||ZB||x1||y1||x2||y2)
 *  S1 = SM3(0x02||yV||inner)
 *  S2 = SM3(0x03||yV||inner)
 * ============================================================ */

static void sm2_kex_compute_inner(uint8_t inner[32],
                                 const sm2_affine_t *V,
                                 const uint8_t ZA[32],
                                 const uint8_t ZB[32],
                                 const sm2_affine_t *RA,
                                 const sm2_affine_t *RB)
{
    uint8_t xV[32], yV[32], x1[32], y1[32], x2[32], y2[32];
    fp_to_bytes(xV, &V->x);
    fp_to_bytes(yV, &V->y);
    fp_to_bytes(x1, &RA->x);
    fp_to_bytes(y1, &RA->y);
    fp_to_bytes(x2, &RB->x);
    fp_to_bytes(y2, &RB->y);

    sm3_ctx c;
    sm3_init(&c);
    sm3_update(&c, xV, 32);
    sm3_update(&c, ZA, 32);
    sm3_update(&c, ZB, 32);
    sm3_update(&c, x1, 32);
    sm3_update(&c, y1, 32);
    sm3_update(&c, x2, 32);
    sm3_update(&c, y2, 32);
    sm3_final(&c, inner);
}

static void sm2_kex_compute_S(uint8_t Sout[32],
                              uint8_t prefix,
                              const sm2_affine_t *V,
                              const uint8_t inner[32])
{
    uint8_t yV[32];
    fp_to_bytes(yV, &V->y);

    sm3_ctx c;
    sm3_init(&c);
    sm3_update(&c, &prefix, 1);
    sm3_update(&c, yV, 32);
    sm3_update(&c, inner, 32);
    sm3_final(&c, Sout);
}

/* ============================================================
 *  协议：发起方 A
 * ============================================================ */

int sm2_kex_initiator_gen_RA(sm2_affine_t *RA,
                             const uint8_t rA[32])
{
    if (!RA || !rA) return 0;

    sm2_affine_t G;
    sm2_get_base_affine(&G);

    sm2_jacobian_t J;
    sm2_scalar_mul_ladder_ct(&J, &G, rA);
    sm2_jacobian_to_affine(RA, &J);

    return sm2_point_is_valid(RA);
}

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
                                  const sm2_affine_t *RB)
{
    if (!K || !S1 || !SA || !dA || !PA || !idA || !rA || !RA || !PB || !idB || !RB) return 0;
    if (!sm2_point_is_valid(PA) || !sm2_point_is_valid(PB) || !sm2_point_is_valid(RA) || !sm2_point_is_valid(RB)) return 0;

    uint8_t ZA[32], ZB[32];
    if (!sm2_kex_compute_Z(ZA, idA, idA_len, PA)) return 0;
    if (!sm2_kex_compute_Z(ZB, idB, idB_len, PB)) return 0;

    /* x1' from RA.x, x2' from RB.x */
    fn_t x1d, x2d;
    sm2_kex_x_dash(&x1d, &RA->x);
    sm2_kex_x_dash(&x2d, &RB->x);

    /* tA = (dA + x1' * rA) mod n */
    uint8_t tA[32];
    if (!sm2_compute_t(tA, dA, rA, &x1d)) return 0;

    /* V = [tA] ( [x2']RB + PB ) */
    sm2_affine_t V;
    if (!sm2_compute_shared_point(&V, tA, &x2d, RB, PB)) return 0;

    /* KDF input: xV||yV||ZA||ZB */
    uint8_t xV[32], yV[32];
    fp_to_bytes(xV, &V.x);
    fp_to_bytes(yV, &V.y);

    uint8_t kdf_in[32+32+32+32];
    memcpy(kdf_in + 0,   xV, 32);
    memcpy(kdf_in + 32,  yV, 32);
    memcpy(kdf_in + 64,  ZA, 32);
    memcpy(kdf_in + 96,  ZB, 32);
    sm3_kdf(K, klen, kdf_in, sizeof(kdf_in));

    /* inner & confirmations */
    uint8_t inner[32];
    sm2_kex_compute_inner(inner, &V, ZA, ZB, RA, RB);

    /* A sends S1 = prefix 0x02 */
    sm2_kex_compute_S(S1, 0x02, &V, inner);
    /* A expects peer to send S2 = prefix 0x03, store locally in SA */
    sm2_kex_compute_S(SA, 0x03, &V, inner);

    if (peer_S2) {
        if (memcmp(peer_S2, SA, 32) != 0) return 0;
    }
    return 1;
}

/* ============================================================
 *  协议：响应方 B
 * ============================================================ */

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
                                  const uint8_t *idA, size_t idA_len)
{
    if (!RB || !K || !S2 || !SB || !peer_S1 || !dB || !PB || !idB || !rB || !RA || !PA || !idA) return 0;
    if (!sm2_point_is_valid(PA) || !sm2_point_is_valid(PB) || !sm2_point_is_valid(RA)) return 0;

    /* generate RB = [rB]G */
    /* 如果调用者已经给了合法 RB，则直接使用；否则内部生成 RB = [rB]G */
    if (!sm2_point_is_valid(RB)) {
        if (!sm2_kex_initiator_gen_RA(RB, rB)) return 0;
    }


    uint8_t ZA[32], ZB[32];
    if (!sm2_kex_compute_Z(ZA, idA, idA_len, PA)) return 0;
    if (!sm2_kex_compute_Z(ZB, idB, idB_len, PB)) return 0;

    /* x2' from RB.x, x1' from RA.x */
    fn_t x1d, x2d;
    sm2_kex_x_dash(&x1d, &RA->x);
    sm2_kex_x_dash(&x2d, &RB->x);

    /* tB = (dB + x2' * rB) mod n */
    uint8_t tB[32];
    if (!sm2_compute_t(tB, dB, rB, &x2d)) return 0;

    /* U = [tB] ( [x1']RA + PA ) */
    sm2_affine_t U;
    if (!sm2_compute_shared_point(&U, tB, &x1d, RA, PA)) return 0;

    /* KDF: xU||yU||ZA||ZB (注意 ZA=发起方, ZB=响应方) */
    uint8_t xU[32], yU[32];
    fp_to_bytes(xU, &U.x);
    fp_to_bytes(yU, &U.y);

    uint8_t kdf_in[32+32+32+32];
    memcpy(kdf_in + 0,   xU, 32);
    memcpy(kdf_in + 32,  yU, 32);
    memcpy(kdf_in + 64,  ZA, 32);
    memcpy(kdf_in + 96,  ZB, 32);
    sm3_kdf(K, klen, kdf_in, sizeof(kdf_in));

    /* inner & confirmations */
    uint8_t inner[32];
    sm2_kex_compute_inner(inner, &U, ZA, ZB, RA, RB);

    /* B verifies peer S1 = prefix 0x02 */
    uint8_t expect_S1[32];
    sm2_kex_compute_S(expect_S1, 0x02, &U, inner);
    if (memcmp(expect_S1, peer_S1, 32) != 0) return 0;

    /* B sends S2 = prefix 0x03 */
    sm2_kex_compute_S(S2, 0x03, &U, inner);

    /* SB：B 侧“本端确认值”（可选留着对照/调试） */
    memcpy(SB, expect_S1, 32);
    return 1;
}


