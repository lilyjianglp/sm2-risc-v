#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "fn.h"
#include "fp.h"
#include "sm2_curve.h"
#include "sm2_scalar.h"
#include "sm2_kex.h"
static inline uint64_t rdcycle(void) {
    uint64_t c;
    asm volatile ("rdcycle %0" : "=r"(c));
    return c;
}

/* ========= 伪随机：仅测试用 ========= */
static uint64_t g_state = 0xC0FFEE123456789AULL;

static uint64_t xorshift64(void) {
    uint64_t x = g_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    g_state = x;
    return x;
}

static void rand_bytes(uint8_t *out, size_t n) {
    for (size_t i = 0; i < n; ) {
        uint64_t r = xorshift64();
        for (int j = 0; j < 8 && i < n; j++, i++) {
            out[i] = (uint8_t)(r >> (56 - 8*j));
        }
    }
}

/* ========= 打印工具 ========= */
static void dump_hex(const char *tag, const uint8_t *p, size_t n) {
    printf("%s (%llu): ", tag, (unsigned long long)n);
    for (size_t i = 0; i < n; i++) printf("%02x", p[i]);
    printf("\n");
}


static void dump_fp(const char *tag, const fp_t *x) {
    uint8_t b[32];
    fp_to_bytes(b, x);
    dump_hex(tag, b, 32);
}

static void dump_point(const char *tag, const sm2_affine_t *P) {
    printf("%s:\n", tag);
    dump_fp("  x", &P->x);
    dump_fp("  y", &P->y);
}

static int bytes_is_zero(const uint8_t a[32]) {
    uint8_t v = 0;
    for (int i = 0; i < 32; i++) v |= a[i];
    return v == 0;
}

/* 生成随机标量 [1, n-1]（用 fn 做 normalize，保证合法） */
static void rand_scalar_mod_n(uint8_t out[32]) {
    uint8_t raw[32];
    fn_t t;

    for (;;) {
        rand_bytes(raw, 32);
        fn_from_be(&t, raw);
        fn_normalize(&t);
        fn_to_be(out, &t);
        if (!bytes_is_zero(out)) break;
    }
}

/* 由私钥 d 计算公钥 P=[d]G（与库内部标量乘保持一致） */
static int derive_public(sm2_affine_t *P, const uint8_t d[32]) {
    sm2_affine_t G;
    sm2_get_base_affine(&G);

    sm2_jacobian_t J;
    sm2_scalar_mul_ladder_ct(&J, &G, d);
    sm2_jacobian_to_affine(P, &J);

    return sm2_point_is_valid(P);
}

 int main(int argc, char **argv) {
    uint64_t start_time = rdcycle();  // 记录程序开始时的时间

    int rounds = 100;
    size_t klen = 32;
    int verbose = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-v")) verbose = 1;
        else if (!strcmp(argv[i], "-n") && i + 1 < argc) rounds = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-k") && i + 1 < argc) klen = (size_t)atoi(argv[++i]);
        else {
            printf("Usage: %s [-n rounds] [-k keylen] [-v]\n", argv[0]);
            return 0;
        }
    }

    printf("[sm2_kex test] rounds=%d, klen=%llu, verbose=%d\n",
       rounds, (unsigned long long)klen, verbose);

    /* 典型 ID（SM2 标准里经常用的占位） */
    const uint8_t idA[] = "Alice";
    const uint8_t idB[] = "Bob";

    for (int tcase = 0; tcase < rounds; tcase++) {
        /* 1) 生成长期密钥对 (dA, PA), (dB, PB) */
        uint8_t dA[32], dB[32];
        rand_scalar_mod_n(dA);
        rand_scalar_mod_n(dB);

        sm2_affine_t PA, PB;
        if (!derive_public(&PA, dA) || !derive_public(&PB, dB)) {
            printf("[FAIL] derive_public at round %d\n", tcase);
            return 1;
        }

        /* 2) 生成临时标量 rA, rB；计算 RA=[rA]G，RB=[rB]G */
        uint8_t rA[32], rB[32];
        rand_scalar_mod_n(rA);
        rand_scalar_mod_n(rB);

        sm2_affine_t RA, RB;
        if (!sm2_kex_initiator_gen_RA(&RA, rA)) {
            printf("[FAIL] gen RA at round %d\n", tcase);
            return 1;
        }

        /* 3) 响应方：计算 Kb、S2，并校验 peer_S1（所以需要先让 A 生成 S1）*/
        uint8_t *KA = (uint8_t*)malloc(klen);
        uint8_t *KB = (uint8_t*)malloc(klen);
        if (!KA || !KB) { printf("[FAIL] malloc\n"); return 1; }

        uint8_t S1[32], SA[32];
        if (!sm2_kex_initiator_gen_RA(&RB, rB)) {
            printf("[FAIL] gen RB at round %d\n", tcase);
            return 1;
        }
        
        if (!sm2_kex_initiator_compute_key(
                KA, klen,
                S1, SA,
                NULL,              /* peer_S2: 暂不校验 */
                dA, &PA, idA, sizeof(idA)-1,
                rA, &RA,
                &PB, idB, sizeof(idB)-1,
                &RB)) {
            printf("[FAIL] initiator_compute_key (phase1) at round %d\n", tcase);
            return 1;
        }

        /* ---- 修正流程：重新跑一遍 ---- */
        if (!sm2_kex_initiator_gen_RA(&RB, rB)) {
            printf("[FAIL] gen RB at round %d\n", tcase);
            return 1;
        }

        if (!sm2_kex_initiator_compute_key(
                KA, klen,
                S1, SA,
                NULL,              /* peer_S2: 暂不校验 */
                dA, &PA, idA, sizeof(idA)-1,
                rA, &RA,
                &PB, idB, sizeof(idB)-1,
                &RB)) {
            printf("[FAIL] initiator_compute_key (phase1) at round %d\n", tcase);
            return 1;
        }

        uint8_t S2[32], SB[32];
        sm2_affine_t RB_out;
        memset(&RB_out, 0, sizeof(RB_out));
        RB_out = RB;

        if (!sm2_kex_responder_compute_key(
                &RB_out,
                KB, klen,
                S2, SB,
                S1,              /* peer_S1 */
                dB, &PB, idB, sizeof(idB)-1,
                rB,
                &RA,
                &PA, idA, sizeof(idA)-1)) {
            printf("[FAIL] responder_compute_key at round %d\n", tcase);
            return 1;
        }

        /* A 侧：用 peer_S2=S2 做最终确认（比较 SA 与 S2） */
        uint8_t S1_2[32], SA_2[32];
        uint8_t *KA2 = (uint8_t*)malloc(klen);
        if (!KA2) { printf("[FAIL] malloc KA2\n"); return 1; }

        if (!sm2_kex_initiator_compute_key(
                KA2, klen,
                S1_2, SA_2,
                S2,              /* peer_S2: 现在校验 */
                dA, &PA, idA, sizeof(idA)-1,
                rA, &RA,
                &PB, idB, sizeof(idB)-1,
                &RB_out)) {
            printf("[FAIL] initiator_compute_key (phase2 verify S2) at round %d\n", tcase);
            return 1;
        }

        /* 4) 检查双方密钥一致 */
        if (memcmp(KA2, KB, klen) != 0) {
            printf("[FAIL] key mismatch at round %d\n", tcase);
            return 1;
        }

        /* 5) 检查 S1/S2 的自洽性：B 侧的 SB 保存的是 expect_S1（你代码里 memcpy(SB, expect_S1)） */
        if (memcmp(S1_2, SB, 32) != 0) {
            printf("[FAIL] S1 mismatch (A vs B expect) at round %d\n", tcase);
            return 1;
        }

        free(KA);
        free(KB);
        free(KA2);
    }

    uint64_t end_time = rdcycle();  // 记录程序结束时的时间
    printf("[CYCLES] Total execution time: %llu cycles\n", (unsigned long long)(end_time - start_time));

    printf("[PASS] sm2_kex end-to-end tests OK ✅\n");
    return 0;
}

