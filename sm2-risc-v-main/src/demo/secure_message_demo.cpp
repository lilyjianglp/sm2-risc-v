#include "sm4_channel.h"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace {

void print_menu()
{
    std::cout << "=========================================\n"
              << " RV-GM-Secure: RISC-V GM Crypto Backend\n"
              << "=========================================\n"
              << "1. Run SM2 key exchange demo\n"
              << "2. Run SM4 secure message demo\n"
              << "3. Run correctness test\n"
              << "4. Run performance benchmark\n"
              << "0. Exit\n"
              << "=========================================\n"
              << "Please select: ";
}

void run_sm2_demo(bool use_optimized_sm2)
{
    std::cout << "\n[SM2 KEX] Running "
              << (use_optimized_sm2 ? "optimized" : "baseline")
              << " SM2 key exchange...\n";

    Sm2KexResult result = run_sm2_key_exchange_demo(use_optimized_sm2);
    bool keys_match = result.success &&
                      !result.shared_key_a.empty() &&
                      result.shared_key_a == result.shared_key_b;

    if (!result.success) {
        std::cout << "[SM2 KEX] Shared key generation failed.\n\n";
        return;
    }

    std::cout << "[SM2 KEX] Shared key generated successfully.\n"
              << "[SM2 KEX] KA == KB: " << (keys_match ? "true" : "false") << "\n"
              << "[SM2 KEX] SharedKey: " << bytes_to_hex(result.shared_key_a) << "\n\n";
}

void run_secure_message_demo(bool use_optimized_sm2)
{
    std::cout << "\nPlease input message:\n";
    std::string message;
    std::getline(std::cin, message);

    std::cout << "\n[SM2 KEX] Running "
              << (use_optimized_sm2 ? "optimized" : "baseline")
              << " SM2 key exchange...\n";

    SecureMessageResult result = run_sm4_secure_message_demo(message, use_optimized_sm2);
    if (!result.kex_success) {
        std::cout << "[SM2 KEX] Shared key generation failed or KA != KB.\n"
                  << "[FAIL] Secure message demo failed.\n\n";
        return;
    }

    std::cout << "[SM2 KEX] Shared key generated successfully.\n"
              << "[SM2 KEX] KA == KB: true\n\n"
              << "[KDF] Deriving SM4 session key from SM2 shared key...\n"
              << "[KDF] SM4 key: " << bytes_to_hex(result.sm4_key_digest) << "\n"
              << "[SM4] IV: " << bytes_to_hex(result.iv) << "\n"
              << "[SM4] Encrypting message...\n"
              << "[SM4] Ciphertext: " << bytes_to_hex(result.ciphertext) << "\n\n"
              << "[SM4] Decrypting message...\n"
              << "[SM4] Plaintext: " << result.decrypted_text << "\n\n"
              << (result.decrypt_success ? "[PASS] Secure message demo success."
                                         : "[FAIL] Secure message demo failed.")
              << "\n\n";
}

} // namespace

int main()
{
#ifdef USE_FP_MONT
    const bool use_optimized_sm2 = true;
#else
    const bool use_optimized_sm2 = false;
#endif

    for (;;) {
        print_menu();
        int choice = -1;
        if (!(std::cin >> choice)) {
            std::cout << "\nInvalid input.\n";
            return 1;
        }
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        switch (choice) {
        case 1:
            run_sm2_demo(use_optimized_sm2);
            break;
        case 2:
            run_secure_message_demo(use_optimized_sm2);
            break;
        case 3:
            std::cout << "\n[INFO] Running correctness test...\n"
                      << "[CMD] make run-pornin-inv-asm USE_FULL_INV_ASM=1\n\n";
            std::system("make run-pornin-inv-asm USE_FULL_INV_ASM=1");
            std::cout << "\n";
            break;
        case 4:
            std::cout << "\n[INFO] Running performance benchmark...\n"
                      << "[CMD] make perf-pornin-asm USE_FULL_INV_ASM=1\n\n";
            std::system("make perf-pornin-asm USE_FULL_INV_ASM=1");
            std::cout << "\n";
            break;
        case 0:
            return 0;
        default:
            std::cout << "\nUnknown selection.\n\n";
            break;
        }
    }
}
