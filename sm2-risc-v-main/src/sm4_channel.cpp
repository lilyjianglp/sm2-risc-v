#include "sm4_channel.h"

extern "C" {
#include "sm2_curve.h"
#include "sm2_kex.h"
#include "sm2_scalar.h"
}

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <ctime>

namespace {

constexpr size_t kSm2SharedKeyLen = 32;
constexpr size_t kSm4BlockLen = 16;
constexpr char kSm4DemoKdfLabel[] = "RV-GM-Secure-SM4-Demo";

void fixed_scalar(uint8_t k[32], uint8_t v)
{
    std::memset(k, v, 32);
    k[31] |= 1;
}

void derive_public_key(sm2_affine_t* pub, const uint8_t d[32])
{
    sm2_affine_t g;
    sm2_get_base_affine(&g);

#ifdef USE_FP_MONT
    sm2_scalar_mul_window_ct_schemeB_mont_to_affine(pub, &g, d);
#else
    sm2_jacobian_t j;
    sm2_scalar_mul_window_ct_schemeB(&j, &g, d);
    sm2_jacobian_to_affine(pub, &j);
#endif
}

Sm2KexResult run_compiled_sm2_kex()
{
    static const uint8_t id_a[] = "Alice";
    static const uint8_t id_b[] = "Bob";

    uint8_t d_a[32], d_b[32], r_a[32], r_b[32];
    fixed_scalar(d_a, 0x11);
    fixed_scalar(d_b, 0x22);
    fixed_scalar(r_a, 0x33);
    fixed_scalar(r_b, 0x44);

    sm2_curve_init_once();
#ifdef USE_FP_MONT
    sm2_curve_init_once_mont();
#endif

    sm2_affine_t p_a, p_b, r_pub_a, r_pub_b;
    derive_public_key(&p_a, d_a);
    derive_public_key(&p_b, d_b);

    uint8_t key_a_first[32], key_a_final[32], key_b[32];
    uint8_t s1_local[32], s1_final[32], s2_local[32];
    uint8_t sa_first[32], sa_final[32], sb[32];

#ifdef USE_FP_MONT
    int ok = sm2_kex_initiator_gen_RA_mont(&r_pub_a, r_a) &&
             sm2_kex_initiator_gen_RA_mont(&r_pub_b, r_b) &&
             sm2_kex_initiator_compute_key_mont(
                 key_a_first, sizeof(key_a_first), s1_local, sa_first, nullptr,
                 d_a, &p_a, id_a, sizeof(id_a) - 1, r_a, &r_pub_a,
                 &p_b, id_b, sizeof(id_b) - 1, &r_pub_b);

    sm2_affine_t r_pub_b_copy = r_pub_b;
    ok = ok && sm2_kex_responder_compute_key_mont(
                   &r_pub_b_copy, key_b, sizeof(key_b), s2_local, sb, s1_local,
                   d_b, &p_b, id_b, sizeof(id_b) - 1, r_b, &r_pub_a,
                   &p_a, id_a, sizeof(id_a) - 1);

    ok = ok && sm2_kex_initiator_compute_key_mont(
                   key_a_final, sizeof(key_a_final), s1_final, sa_final,
                   s2_local, d_a, &p_a, id_a, sizeof(id_a) - 1, r_a,
                   &r_pub_a, &p_b, id_b, sizeof(id_b) - 1, &r_pub_b);
#else
    int ok = sm2_kex_initiator_gen_RA(&r_pub_a, r_a) &&
             sm2_kex_initiator_gen_RA(&r_pub_b, r_b) &&
             sm2_kex_initiator_compute_key(
                 key_a_first, sizeof(key_a_first), s1_local, sa_first, nullptr,
                 d_a, &p_a, id_a, sizeof(id_a) - 1, r_a, &r_pub_a,
                 &p_b, id_b, sizeof(id_b) - 1, &r_pub_b);

    sm2_affine_t r_pub_b_copy = r_pub_b;
    ok = ok && sm2_kex_responder_compute_key(
                   &r_pub_b_copy, key_b, sizeof(key_b), s2_local, sb, s1_local,
                   d_b, &p_b, id_b, sizeof(id_b) - 1, r_b, &r_pub_a,
                   &p_a, id_a, sizeof(id_a) - 1);

    ok = ok && sm2_kex_initiator_compute_key(
                   key_a_final, sizeof(key_a_final), s1_final, sa_final,
                   s2_local, d_a, &p_a, id_a, sizeof(id_a) - 1, r_a,
                   &r_pub_a, &p_b, id_b, sizeof(id_b) - 1, &r_pub_b);
#endif

    ok = ok &&
         std::memcmp(key_a_first, key_b, sizeof(key_b)) == 0 &&
         std::memcmp(key_a_final, key_b, sizeof(key_b)) == 0 &&
         std::memcmp(sa_first, s2_local, sizeof(s2_local)) == 0 &&
         std::memcmp(s1_final, s1_local, sizeof(s1_local)) == 0;

    Sm2KexResult result;
    result.success = ok != 0;
    if (result.success) {
        result.shared_key_a.assign(key_a_final, key_a_final + sizeof(key_a_final));
        result.shared_key_b.assign(key_b, key_b + sizeof(key_b));
    }
    return result;
}

uint32_t rotl32(uint32_t x, uint32_t n)
{
    n &= 31;
    return (x << n) | (x >> ((32 - n) & 31));
}

uint32_t load_be32(const uint8_t* p)
{
    return (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) |
           (uint32_t(p[2]) << 8) | uint32_t(p[3]);
}

void store_be32(uint8_t* p, uint32_t x)
{
    p[0] = uint8_t(x >> 24);
    p[1] = uint8_t(x >> 16);
    p[2] = uint8_t(x >> 8);
    p[3] = uint8_t(x);
}

uint32_t sm3_p0(uint32_t x) { return x ^ rotl32(x, 9) ^ rotl32(x, 17); }
uint32_t sm3_p1(uint32_t x) { return x ^ rotl32(x, 15) ^ rotl32(x, 23); }

uint32_t sm3_ff(uint32_t x, uint32_t y, uint32_t z, int j)
{
    return j < 16 ? (x ^ y ^ z) : ((x & y) | (x & z) | (y & z));
}

uint32_t sm3_gg(uint32_t x, uint32_t y, uint32_t z, int j)
{
    return j < 16 ? (x ^ y ^ z) : ((x & y) | ((~x) & z));
}

class Sm3 {
public:
    Sm3()
    {
        h_ = {0x7380166Fu, 0x4914B2B9u, 0x172442D7u, 0xDA8A0600u,
              0xA96F30BCu, 0x163138AAu, 0xE38DEE4Du, 0xB0FB0E4Eu};
    }

    void update(const uint8_t* in, size_t len)
    {
        nbits_ += uint64_t(len) * 8;
        while (len > 0) {
            size_t n = std::min(len, block_.size() - block_len_);
            std::memcpy(block_.data() + block_len_, in, n);
            block_len_ += n;
            in += n;
            len -= n;
            if (block_len_ == block_.size()) {
                compress(block_.data());
                block_len_ = 0;
            }
        }
    }

    std::array<uint8_t, 32> final()
    {
        uint8_t pad[72] = {0x80};
        size_t zero_len = block_len_ < 56 ? (56 - block_len_ - 1)
                                          : (64 + 56 - block_len_ - 1);
        uint64_t bits = nbits_;
        store_be32(pad + 1 + zero_len, uint32_t(bits >> 32));
        store_be32(pad + 1 + zero_len + 4, uint32_t(bits));
        update(pad, 1 + zero_len + 8);

        std::array<uint8_t, 32> out{};
        for (size_t i = 0; i < h_.size(); ++i) {
            store_be32(out.data() + i * 4, h_[i]);
        }
        return out;
    }

private:
    void compress(const uint8_t block[64])
    {
        uint32_t w[68], w1[64];
        for (int i = 0; i < 16; ++i) {
            w[i] = load_be32(block + 4 * i);
        }
        for (int j = 16; j < 68; ++j) {
            uint32_t x = w[j - 16] ^ w[j - 9] ^ rotl32(w[j - 3], 15);
            w[j] = sm3_p1(x) ^ rotl32(w[j - 13], 7) ^ w[j - 6];
        }
        for (int j = 0; j < 64; ++j) {
            w1[j] = w[j] ^ w[j + 4];
        }

        uint32_t a = h_[0], b = h_[1], c = h_[2], d = h_[3];
        uint32_t e = h_[4], f = h_[5], g = h_[6], hh = h_[7];
        for (int j = 0; j < 64; ++j) {
            uint32_t tj = j < 16 ? 0x79CC4519u : 0x7A879D8Au;
            uint32_t ss1 = rotl32(rotl32(a, 12) + e + rotl32(tj, uint32_t(j)), 7);
            uint32_t ss2 = ss1 ^ rotl32(a, 12);
            uint32_t tt1 = sm3_ff(a, b, c, j) + d + ss2 + w1[j];
            uint32_t tt2 = sm3_gg(e, f, g, j) + hh + ss1 + w[j];
            d = c;
            c = rotl32(b, 9);
            b = a;
            a = tt1;
            hh = g;
            g = rotl32(f, 19);
            f = e;
            e = sm3_p0(tt2);
        }

        h_[0] ^= a; h_[1] ^= b; h_[2] ^= c; h_[3] ^= d;
        h_[4] ^= e; h_[5] ^= f; h_[6] ^= g; h_[7] ^= hh;
    }

    std::array<uint32_t, 8> h_{};
    uint64_t nbits_ = 0;
    std::array<uint8_t, 64> block_{};
    size_t block_len_ = 0;
};

std::array<uint8_t, 32> sm3_hash(const uint8_t* data, size_t len)
{
    Sm3 sm3;
    sm3.update(data, len);
    return sm3.final();
}

std::vector<uint8_t> pkcs7_pad(const std::string& plaintext)
{
    std::vector<uint8_t> out(plaintext.begin(), plaintext.end());
    size_t pad_len = kSm4BlockLen - (out.size() % kSm4BlockLen);
    if (pad_len == 0) {
        pad_len = kSm4BlockLen;
    }
    out.insert(out.end(), pad_len, uint8_t(pad_len));
    return out;
}

std::string pkcs7_unpad(const std::vector<uint8_t>& padded)
{
    if (padded.empty() || padded.size() % kSm4BlockLen != 0) {
        throw std::runtime_error("invalid padded plaintext size");
    }
    uint8_t pad_len = padded.back();
    if (pad_len == 0 || pad_len > kSm4BlockLen || pad_len > padded.size()) {
        throw std::runtime_error("invalid PKCS#7 padding");
    }
    for (size_t i = padded.size() - pad_len; i < padded.size(); ++i) {
        if (padded[i] != pad_len) {
            throw std::runtime_error("invalid PKCS#7 padding bytes");
        }
    }
    return std::string(padded.begin(), padded.end() - pad_len);
}

constexpr uint8_t kSm4Sbox[256] = {
    0xd6,0x90,0xe9,0xfe,0xcc,0xe1,0x3d,0xb7,0x16,0xb6,0x14,0xc2,0x28,0xfb,0x2c,0x05,
    0x2b,0x67,0x9a,0x76,0x2a,0xbe,0x04,0xc3,0xaa,0x44,0x13,0x26,0x49,0x86,0x06,0x99,
    0x9c,0x42,0x50,0xf4,0x91,0xef,0x98,0x7a,0x33,0x54,0x0b,0x43,0xed,0xcf,0xac,0x62,
    0xe4,0xb3,0x1c,0xa9,0xc9,0x08,0xe8,0x95,0x80,0xdf,0x94,0xfa,0x75,0x8f,0x3f,0xa6,
    0x47,0x07,0xa7,0xfc,0xf3,0x73,0x17,0xba,0x83,0x59,0x3c,0x19,0xe6,0x85,0x4f,0xa8,
    0x68,0x6b,0x81,0xb2,0x71,0x64,0xda,0x8b,0xf8,0xeb,0x0f,0x4b,0x70,0x56,0x9d,0x35,
    0x1e,0x24,0x0e,0x5e,0x63,0x58,0xd1,0xa2,0x25,0x22,0x7c,0x3b,0x01,0x21,0x78,0x87,
    0xd4,0x00,0x46,0x57,0x9f,0xd3,0x27,0x52,0x4c,0x36,0x02,0xe7,0xa0,0xc4,0xc8,0x9e,
    0xea,0xbf,0x8a,0xd2,0x40,0xc7,0x38,0xb5,0xa3,0xf7,0xf2,0xce,0xf9,0x61,0x15,0xa1,
    0xe0,0xae,0x5d,0xa4,0x9b,0x34,0x1a,0x55,0xad,0x93,0x32,0x30,0xf5,0x8c,0xb1,0xe3,
    0x1d,0xf6,0xe2,0x2e,0x82,0x66,0xca,0x60,0xc0,0x29,0x23,0xab,0x0d,0x53,0x4e,0x6f,
    0xd5,0xdb,0x37,0x45,0xde,0xfd,0x8e,0x2f,0x03,0xff,0x6a,0x72,0x6d,0x6c,0x5b,0x51,
    0x8d,0x1b,0xaf,0x92,0xbb,0xdd,0xbc,0x7f,0x11,0xd9,0x5c,0x41,0x1f,0x10,0x5a,0xd8,
    0x0a,0xc1,0x31,0x88,0xa5,0xcd,0x7b,0xbd,0x2d,0x74,0xd0,0x12,0xb8,0xe5,0xb4,0xb0,
    0x89,0x69,0x97,0x4a,0x0c,0x96,0x77,0x7e,0x65,0xb9,0xf1,0x09,0xc5,0x6e,0xc6,0x84,
    0x18,0xf0,0x7d,0xec,0x3a,0xdc,0x4d,0x20,0x79,0xee,0x5f,0x3e,0xd7,0xcb,0x39,0x48
};

constexpr uint32_t kSm4Fk[4] = {0xa3b1bac6u, 0x56aa3350u, 0x677d9197u, 0xb27022dcu};

constexpr uint32_t kSm4Ck[32] = {
    0x00070e15u,0x1c232a31u,0x383f464du,0x545b6269u,
    0x70777e85u,0x8c939aa1u,0xa8afb6bdu,0xc4cbd2d9u,
    0xe0e7eef5u,0xfc030a11u,0x181f262du,0x343b4249u,
    0x50575e65u,0x6c737a81u,0x888f969du,0xa4abb2b9u,
    0xc0c7ced5u,0xdce3eaf1u,0xf8ff060du,0x141b2229u,
    0x30373e45u,0x4c535a61u,0x686f767du,0x848b9299u,
    0xa0a7aeb5u,0xbcc3cad1u,0xd8dfe6edu,0xf4fb0209u,
    0x10171e25u,0x2c333a41u,0x484f565du,0x646b7279u
};

uint32_t sm4_tau(uint32_t x)
{
    uint8_t b0 = kSm4Sbox[(x >> 24) & 0xff];
    uint8_t b1 = kSm4Sbox[(x >> 16) & 0xff];
    uint8_t b2 = kSm4Sbox[(x >> 8) & 0xff];
    uint8_t b3 = kSm4Sbox[x & 0xff];
    return (uint32_t(b0) << 24) | (uint32_t(b1) << 16) |
           (uint32_t(b2) << 8) | uint32_t(b3);
}

uint32_t sm4_l(uint32_t x)
{
    return x ^ rotl32(x, 2) ^ rotl32(x, 10) ^ rotl32(x, 18) ^ rotl32(x, 24);
}

uint32_t sm4_l_key(uint32_t x)
{
    return x ^ rotl32(x, 13) ^ rotl32(x, 23);
}

std::array<uint32_t, 32> sm4_round_keys(const std::vector<uint8_t>& key)
{
    if (key.size() != kSm4BlockLen) {
        throw std::runtime_error("SM4 key must be 16 bytes");
    }

    uint32_t k[36];
    for (int i = 0; i < 4; ++i) {
        k[i] = load_be32(key.data() + 4 * i) ^ kSm4Fk[i];
    }

    std::array<uint32_t, 32> rk{};
    for (int i = 0; i < 32; ++i) {
        k[i + 4] = k[i] ^ sm4_l_key(sm4_tau(k[i + 1] ^ k[i + 2] ^ k[i + 3] ^ kSm4Ck[i]));
        rk[i] = k[i + 4];
    }
    return rk;
}

void sm4_crypt_block(uint8_t out[16], const uint8_t in[16], const std::array<uint32_t, 32>& rk)
{
    uint32_t x[36];
    for (int i = 0; i < 4; ++i) {
        x[i] = load_be32(in + 4 * i);
    }
    for (int i = 0; i < 32; ++i) {
        x[i + 4] = x[i] ^ sm4_l(sm4_tau(x[i + 1] ^ x[i + 2] ^ x[i + 3] ^ rk[i]));
    }
    store_be32(out, x[35]);
    store_be32(out + 4, x[34]);
    store_be32(out + 8, x[33]);
    store_be32(out + 12, x[32]);
}

} // namespace

Sm2KexResult run_sm2_key_exchange_demo(bool use_optimized_sm2)
{
    (void)use_optimized_sm2;
    return run_compiled_sm2_kex();
}

std::vector<uint8_t> derive_sm4_key(const std::vector<uint8_t>& shared_key)
{
    std::vector<uint8_t> input(shared_key);
    input.insert(input.end(), kSm4DemoKdfLabel,
                 kSm4DemoKdfLabel + std::strlen(kSm4DemoKdfLabel));

    auto digest = sm3_hash(input.data(), input.size());
    return std::vector<uint8_t>(digest.begin(), digest.begin() + kSm4BlockLen);
}

std::vector<uint8_t> generate_random_iv()
{
    std::vector<uint8_t> iv(kSm4BlockLen);

    FILE* urandom = std::fopen("/dev/urandom", "rb");
    if (urandom != nullptr) {
        size_t nread = std::fread(iv.data(), 1, iv.size(), urandom);
        std::fclose(urandom);
        if (nread == iv.size()) {
            return iv;
        }
    }

    uint64_t state = uint64_t(std::time(nullptr)) ^
                     (uint64_t(reinterpret_cast<uintptr_t>(iv.data())) << 1);
    for (auto& byte : iv) {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        byte = uint8_t(state >> 56);
    }
    return iv;
}

std::vector<uint8_t> sm4_encrypt_message(
    const std::vector<uint8_t>& key,
    const std::vector<uint8_t>& iv,
    const std::string& plaintext)
{
    if (iv.size() != kSm4BlockLen) {
        throw std::runtime_error("SM4-CBC IV must be 16 bytes");
    }

    auto rk = sm4_round_keys(key);
    auto padded = pkcs7_pad(plaintext);
    std::vector<uint8_t> ciphertext(padded.size());
    std::array<uint8_t, 16> prev{};
    std::copy(iv.begin(), iv.end(), prev.begin());

    for (size_t off = 0; off < padded.size(); off += kSm4BlockLen) {
        uint8_t block[16];
        for (size_t i = 0; i < kSm4BlockLen; ++i) {
            block[i] = padded[off + i] ^ prev[i];
        }
        sm4_crypt_block(ciphertext.data() + off, block, rk);
        std::copy(ciphertext.begin() + off, ciphertext.begin() + off + kSm4BlockLen, prev.begin());
    }

    return ciphertext;
}

std::string sm4_decrypt_message(
    const std::vector<uint8_t>& key,
    const std::vector<uint8_t>& iv,
    const std::vector<uint8_t>& ciphertext)
{
    if (iv.size() != kSm4BlockLen) {
        throw std::runtime_error("SM4-CBC IV must be 16 bytes");
    }
    if (ciphertext.empty() || ciphertext.size() % kSm4BlockLen != 0) {
        throw std::runtime_error("SM4-CBC ciphertext must be a non-empty block multiple");
    }

    auto rk = sm4_round_keys(key);
    std::reverse(rk.begin(), rk.end());

    std::vector<uint8_t> padded(ciphertext.size());
    std::array<uint8_t, 16> prev{};
    std::copy(iv.begin(), iv.end(), prev.begin());

    for (size_t off = 0; off < ciphertext.size(); off += kSm4BlockLen) {
        uint8_t block[16];
        sm4_crypt_block(block, ciphertext.data() + off, rk);
        for (size_t i = 0; i < kSm4BlockLen; ++i) {
            padded[off + i] = block[i] ^ prev[i];
        }
        std::copy(ciphertext.begin() + off, ciphertext.begin() + off + kSm4BlockLen, prev.begin());
    }

    return pkcs7_unpad(padded);
}

SecureMessageResult run_sm4_secure_message_demo(
    const std::string& message,
    bool use_optimized_sm2)
{
    SecureMessageResult result;
    result.plaintext = message;

    Sm2KexResult kex = run_sm2_key_exchange_demo(use_optimized_sm2);
    result.kex_success = kex.success &&
                         kex.shared_key_a.size() == kSm2SharedKeyLen &&
                         kex.shared_key_a == kex.shared_key_b;
    if (!result.kex_success) {
        return result;
    }

    result.shared_key_digest = kex.shared_key_a;
    result.sm4_key_digest = derive_sm4_key(kex.shared_key_a);
    result.iv = generate_random_iv();
    result.ciphertext = sm4_encrypt_message(result.sm4_key_digest, result.iv, message);
    result.decrypted_text = sm4_decrypt_message(result.sm4_key_digest, result.iv, result.ciphertext);
    result.decrypt_success = result.decrypted_text == message;
    return result;
}

std::string bytes_to_hex(const std::vector<uint8_t>& bytes)
{
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (uint8_t byte : bytes) {
        oss << std::setw(2) << int(byte);
    }
    return oss.str();
}
