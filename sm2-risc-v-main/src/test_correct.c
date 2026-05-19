/*
 * test_correct.c
 * ------------------------------------------------------------
 * SM2 correctness test for current build variant.
 *
 * This file supports all four build variants:
 *   baseline-c    : Prime route + fp C + fn C
 *   baseline-asm  : Prime route + fp ASM + fn ASM
 *   mont-c        : Montgomery route + fp C + fn C
 *   mont-asm      : Montgomery route + fp ASM + fn ASM
 *
 * It does NOT link normal and Montgomery implementations together.
 * Instead, it tests the currently compiled implementation.
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stddef.h>

#include "fn.h"
#include "fp.h"
#include "sm2_curve.h"
#include "sm2_scalar.h"
#include "sm2_kex.h"

/* ============================================================
 * Helpers
 * ============================================================ */

static int bytes_equal(const uint8_t *a, const uint8_t *b, size_t len) {
    return memcmp(a, b, len) == 0;
}

static void print_hex(const char *label, const uint8_t *data, size_t len) {
    printf("%s = ", label);
    for (size_t i = 0; i < len; i++) {
        printf("%02x", data[i]);
    }
    printf("\n");
}

static void fixed_scalar(uint8_t k[32], uint8_t v) {
    memset(k, v, 32);
    k[31] |= 1;
}

static int is_zero_scalar(const uint8_t k[32]) {
    uint8_t x = 0;
    for (int i = 0; i < 32; i++) {
        x |= k[i];
    }
    return x == 0;
}

/* deterministic RNG for reproducible random tests */
static uint64_t g_test_rng_state = 0x123456789abcdef0ULL;

static uint64_t test_xorshift64(void) {
    uint64_t x = g_test_rng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    g_test_rng_state = x;
    return x;
}

static void scalar_random_nonzero_test(uint8_t k[32]) {
    do {
        for (int i = 0; i < 4; i++) {
            uint64_t w = test_xorshift64();
            for (int j = 0; j < 8; j++) {
                k[i * 8 + j] = (uint8_t)(w >> (56 - 8 * j));
            }
        }

        /*
         * Keep the top bit clear for stable test scalars.
         * This is not cryptographic randomness; it is only for tests.
         */
        k[0] &= 0x7f;
    } while (is_zero_scalar(k));
}

/* ============================================================
 * Build variant wrappers
 * ============================================================ */

static const char *build_variant_name(void) {
#ifdef USE_FP_MONT
#ifdef USE_RV64_ASM_MUL
    return "mont-asm";
#else
    return "mont-c";
#endif
#else
#ifdef USE_ASM_FP
    return "baseline-asm";
#else
    return "baseline-c";
#endif
#endif
}

static void derive_public_key(sm2_affine_t *Ppub, const uint8_t d[32]) {
    sm2_affine_t G;
    sm2_get_base_affine(&G);

#ifdef USE_FP_MONT
    sm2_scalar_mul_window_ct_schemeB_mont_to_affine(Ppub, &G, d);
#else
    sm2_jacobian_t J;
    sm2_scalar_mul_window_ct_schemeB(&J, &G, d);
    sm2_jacobian_to_affine(Ppub, &J);
#endif
}

static int kex_gen_R(sm2_affine_t *R, const uint8_t r[32]) {
#ifdef USE_FP_MONT
    return sm2_kex_initiator_gen_RA_mont(R, r);
#else
    return sm2_kex_initiator_gen_RA(R, r);
#endif
}

static int kex_initiator_compute(uint8_t K[32],
                                 uint8_t S1[32],
                                 uint8_t SA[32],
                                 const uint8_t *S2_in,
                                 const uint8_t dA[32],
                                 const sm2_affine_t *PA,
                                 const uint8_t *idA,
                                 size_t idA_len,
                                 const uint8_t rA[32],
                                 const sm2_affine_t *RA,
                                 const sm2_affine_t *PB,
                                 const uint8_t *idB,
                                 size_t idB_len,
                                 const sm2_affine_t *RB) {
#ifdef USE_FP_MONT
    return sm2_kex_initiator_compute_key_mont(
        K, 32, S1, SA, S2_in,
        dA, PA, idA, idA_len,
        rA, RA,
        PB, idB, idB_len,
        RB
    );
#else
    return sm2_kex_initiator_compute_key(
        K, 32, S1, SA, S2_in,
        dA, PA, idA, idA_len,
        rA, RA,
        PB, idB, idB_len,
        RB
    );
#endif
}

static int kex_responder_compute(sm2_affine_t *RB,
                                 uint8_t K[32],
                                 uint8_t S2[32],
                                 uint8_t SB[32],
                                 const uint8_t S1_in[32],
                                 const uint8_t dB[32],
                                 const sm2_affine_t *PB,
                                 const uint8_t *idB,
                                 size_t idB_len,
                                 const uint8_t rB[32],
                                 const sm2_affine_t *RA,
                                 const sm2_affine_t *PA,
                                 const uint8_t *idA,
                                 size_t idA_len) {
#ifdef USE_FP_MONT
    return sm2_kex_responder_compute_key_mont(
        RB, K, 32, S2, SB, S1_in,
        dB, PB, idB, idB_len,
        rB, RA, PA, idA, idA_len
    );
#else
    return sm2_kex_responder_compute_key(
        RB, K, 32, S2, SB, S1_in,
        dB, PB, idB, idB_len,
        rB, RA, PA, idA, idA_len
    );
#endif
}

/* ============================================================
 * Single SM2 key exchange correctness case
 * ============================================================ */

static int run_one_kex_case(const char *name,
                            const uint8_t dA[32],
                            const uint8_t dB[32],
                            const uint8_t rA[32],
                            const uint8_t rB[32],
                            int print_result) {
    static const uint8_t idA[] = "Alice";
    static const uint8_t idB[] = "Bob";

    sm2_affine_t PA, PB;
    sm2_affine_t RA, RB;

    uint8_t KA1[32], KA2[32], KB[32];
    uint8_t S1_first[32], S1_final[32], S2_resp[32];
    uint8_t SA1[32], SA2[32], SB[32];

    derive_public_key(&PA, dA);
    derive_public_key(&PB, dB);

    if (!kex_gen_R(&RA, rA)) {
        printf("FAIL: %s: generate RA\n", name);
        return 0;
    }

    if (!kex_gen_R(&RB, rB)) {
        printf("FAIL: %s: generate RB\n", name);
        return 0;
    }

    /*
     * Initiator first pass:
     * produces KA1 and S1_first for responder verification.
     */
    if (!kex_initiator_compute(
            KA1, S1_first, SA1, NULL,
            dA, &PA, idA, sizeof(idA) - 1,
            rA, &RA,
            &PB, idB, sizeof(idB) - 1,
            &RB)) {
        printf("FAIL: %s: initiator first pass\n", name);
        return 0;
    }

    /*
     * Responder pass:
     * verifies S1_first and produces KB and S2_resp.
     */
    {
        sm2_affine_t RB_copy = RB;

        if (!kex_responder_compute(
                &RB_copy, KB, S2_resp, SB, S1_first,
                dB, &PB, idB, sizeof(idB) - 1,
                rB, &RA, &PA, idA, sizeof(idA) - 1)) {
            printf("FAIL: %s: responder compute\n", name);
            return 0;
        }
    }

    /*
     * Initiator second pass:
     * verifies S2_resp and produces final KA2.
     */
    if (!kex_initiator_compute(
            KA2, S1_final, SA2, S2_resp,
            dA, &PA, idA, sizeof(idA) - 1,
            rA, &RA,
            &PB, idB, sizeof(idB) - 1,
            &RB)) {
        printf("FAIL: %s: initiator verify pass\n", name);
        return 0;
    }

    if (!bytes_equal(KA1, KB, 32)) {
        printf("FAIL: %s: KA1 != KB\n", name);
        print_hex("KA1", KA1, 32);
        print_hex("KB ", KB, 32);
        return 0;
    }

    if (!bytes_equal(KA2, KB, 32)) {
        printf("FAIL: %s: KA2 != KB\n", name);
        print_hex("KA2", KA2, 32);
        print_hex("KB ", KB, 32);
        return 0;
    }

    if (!bytes_equal(S1_final, S1_first, 32)) {
        printf("FAIL: %s: S1_final != S1_first\n", name);
        print_hex("S1_first", S1_first, 32);
        print_hex("S1_final", S1_final, 32);
        return 0;
    }

    if (print_result) {
        printf("PASS: %s\n", name);
        print_hex("SharedKey", KB, 32);
        print_hex("S1", S1_final, 32);
        print_hex("S2", S2_resp, 32);
    } else {
        printf("PASS: %s\n", name);
    }

    return 1;
}

/* ============================================================
 * Correctness suite
 * ============================================================ */

static int run_fixed_kex_test(void) {
    static const uint8_t dA[32] = {
        0x12,0x34,0x56,0x78,0x9a,0xbc,0xde,0xf0,
        0x11,0x22,0x33,0x44,0x55,0x66,0x77,0x88,
        0x99,0xaa,0xbb,0xcc,0xdd,0xee,0xf1,0x02,
        0x13,0x24,0x35,0x46,0x57,0x68,0x79,0x8a
    };

    static const uint8_t dB[32] = {
        0x23,0x45,0x67,0x89,0xab,0xcd,0xef,0x01,
        0x10,0x32,0x54,0x76,0x98,0xba,0xdc,0xfe,
        0x12,0x34,0x56,0x78,0x87,0x65,0x43,0x21,
        0x0f,0x1e,0x2d,0x3c,0x4b,0x5a,0x69,0x78
    };

    static const uint8_t rA[32] = {
        0x34,0x56,0x78,0x9a,0xbc,0xde,0xf0,0x12,
        0x21,0x43,0x65,0x87,0xa9,0xcb,0xed,0x0f,
        0x13,0x57,0x9b,0xdf,0x24,0x68,0xac,0xe0,
        0x10,0x20,0x30,0x40,0x50,0x60,0x70,0x80
    };

    static const uint8_t rB[32] = {
        0x45,0x67,0x89,0xab,0xcd,0xef,0x01,0x23,
        0x32,0x54,0x76,0x98,0xba,0xdc,0xfe,0x10,
        0x24,0x68,0xac,0xe0,0x35,0x79,0xbd,0xf1,
        0x11,0x22,0x33,0x44,0x55,0x66,0x77,0x88
    };

    return run_one_kex_case("fixed SM2 KEX test", dA, dB, rA, rB, 1);
}

static int run_random_kex_tests(int rounds) {
    for (int i = 0; i < rounds; i++) {
        uint8_t dA[32], dB[32], rA[32], rB[32];
        char name[64];

        scalar_random_nonzero_test(dA);
        scalar_random_nonzero_test(dB);
        scalar_random_nonzero_test(rA);
        scalar_random_nonzero_test(rB);

        snprintf(name, sizeof(name), "random SM2 KEX test #%d", i + 1);

        if (!run_one_kex_case(name, dA, dB, rA, rB, 0)) {
            return 0;
        }
    }

    printf("PASS: random SM2 KEX tests (%d rounds)\n", rounds);
    return 1;
}

static int run_correctness_suite(void) {
    printf("\n=== SM2 correctness suite (%s) ===\n", build_variant_name());

    if (!run_fixed_kex_test()) {
        return 0;
    }

    if (!run_random_kex_tests(10)) {
        return 0;
    }

    printf("=== SM2 correctness suite PASSED (%s) ===\n\n", build_variant_name());
    return 1;
}

/* ============================================================
 * main
 * ============================================================ */

int main(void) {
    printf("SM2 Correctness Test\n");
    printf("====================\n");
    printf("Build variant: %s\n", build_variant_name());

    sm2_curve_init_once();

    if (!run_correctness_suite()) {
        puts("[FAIL] SM2 correctness test failed");
        return 1;
    }

    puts("[PASS] SM2 correctness test passed");
    return 0;
}