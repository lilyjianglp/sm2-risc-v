#pragma once

#include <cstdint>
#include <string>
#include <vector>

struct Sm2KexResult {
    bool success = false;
    std::vector<uint8_t> shared_key_a;
    std::vector<uint8_t> shared_key_b;
};

struct SecureMessageResult {
    bool kex_success = false;
    bool decrypt_success = false;

    std::string plaintext;
    std::string decrypted_text;

    std::vector<uint8_t> shared_key_digest;
    std::vector<uint8_t> sm4_key_digest;
    std::vector<uint8_t> iv;
    std::vector<uint8_t> ciphertext;
};

Sm2KexResult run_sm2_key_exchange_demo(bool use_optimized_sm2);

std::vector<uint8_t> derive_sm4_key(
    const std::vector<uint8_t>& shared_key
);

std::vector<uint8_t> generate_random_iv();

std::vector<uint8_t> sm4_encrypt_message(
    const std::vector<uint8_t>& key,
    const std::vector<uint8_t>& iv,
    const std::string& plaintext
);

std::string sm4_decrypt_message(
    const std::vector<uint8_t>& key,
    const std::vector<uint8_t>& iv,
    const std::vector<uint8_t>& ciphertext
);

SecureMessageResult run_sm4_secure_message_demo(
    const std::string& message,
    bool use_optimized_sm2
);

std::string bytes_to_hex(const std::vector<uint8_t>& bytes);
