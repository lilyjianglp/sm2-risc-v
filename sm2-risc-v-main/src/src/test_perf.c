/*
 * test_perf.c
 * ---------------------------------------
 * SM2 perf profiling workload (no rdcycle)
 *
 * Purpose:
 *   - Provide clean workloads for Linux perf stat / perf record
 *   - Avoid mixing perf sampling with in-program rdcycle microbenchmarking
 *
 * Modes:
 *   genra       : profile sm2_kex_initiator_gen_RA_mont()
 *   scalar      : profile fixed-base scalar multiplication core
 *   responder   : profile sm2_kex_responder_compute_key_mont()
 *   initiator   : profile sm2_kex_initiator_compute_key_mont() verify path
 *   full        : profile responder + initiator path in one loop
 *   all         : run all above workloads sequentially
 *
 * Usage examples:
 *   ./test_perf
 *   ./test_perf responder
 *   ./test_perf scalar 50000
 *   perf stat ./test_perf responder 20000
 *   perf record ./test_perf scalar 50000
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "fn.h"
#include "fp.h"
#include "sm2_curve.h"
#include "sm2_scalar.h"
#include "sm2_kex.h"

#ifndef USE_FP_MONT
# error "test_perf.c is intended for Montgomery build. Compile with USE_FP_MONT enabled."
#endif

#ifndef DEFAULT_PERF_ITERS
#define DEFAULT_PERF_ITERS 20000UL
#endif

volatile uint8_t bench_sink = 0;

static inline void sink_mem(const void *p, size_t n) {
    const volatile uint8_t *b = (const volatile uint8_t *)p;
    for (size_t i = 0; i < n; i++) bench_sink ^= b[i];
}

static void fixed_scalar(uint8_t k[32], uint8_t v) {
    memset(k, v, 32);
    k[31] |= 1;
}

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [mode] [iters]\n"
        "\n"
        "Modes:\n"
        "  genra       Profile sm2_kex_initiator_gen_RA_mont()\n"
        "  scalar      Profile fixed-base scalar multiplication core\n"
        "  responder   Profile sm2_kex_responder_compute_key_mont()\n"
        "  initiator   Profile sm2_kex_initiator_compute_key_mont()\n"
        "  full        Profile responder + initiator sequence\n"
        "  all         Run all workloads sequentially (default)\n"
        "\n"
        "Examples:\n"
        "  %s responder\n"
        "  %s scalar 50000\n",
        prog, prog, prog);
}

typedef struct {
    sm2_affine_t G;
    sm2_affine_t PA, PB;
    sm2_affine_t RA, RB;
    uint8_t dA[32], dB[32], rA[32], rB[32];
    uint8_t KA[32], KB[32];
    uint8_t S1[32], S2[32], SA[32], SB[32];
    const uint8_t *idA;
    size_t idA_len;
    const uint8_t *idB;
    size_t idB_len;
} perf_ctx_t;

static int setup_perf_ctx(perf_ctx_t *ctx) {
    static const uint8_t idA_local[] = "Alice";
    static const uint8_t idB_local[] = "Bob";

    memset(ctx, 0, sizeof(*ctx));

    ctx->idA = idA_local;
    ctx->idA_len = sizeof(idA_local) - 1;
    ctx->idB = idB_local;
    ctx->idB_len = sizeof(idB_local) - 1;

    sm2_curve_init_once();
    sm2_get_base_affine(&ctx->G);

    fixed_scalar(ctx->dA, 0x11);
    fixed_scalar(ctx->dB, 0x22);
    fixed_scalar(ctx->rA, 0x33);
    fixed_scalar(ctx->rB, 0x44);

    sm2_scalar_mul_window_ct_schemeB_mont_to_affine(&ctx->PA, &ctx->G, ctx->dA);
    sm2_scalar_mul_window_ct_schemeB_mont_to_affine(&ctx->PB, &ctx->G, ctx->dB);
    if (!sm2_kex_initiator_gen_RA_mont(&ctx->RA, ctx->rA)) {
        fprintf(stderr, "setup failed: gen RA\n");
        return 0;
    }
    if (!sm2_kex_initiator_gen_RA_mont(&ctx->RB, ctx->rB)) {
        fprintf(stderr, "setup failed: gen RB\n");
        return 0;
    }

    /*
     * Prepare valid protocol-side inputs:
     *   1) initiator first pass -> produces S1 for responder input
     *   2) responder pass       -> produces S2 for initiator verify input
     */
    if (!sm2_kex_initiator_compute_key_mont(
            ctx->KA, 32, ctx->S1, ctx->SA, NULL,
            ctx->dA, &ctx->PA, ctx->idA, ctx->idA_len,
            ctx->rA, &ctx->RA,
            &ctx->PB, ctx->idB, ctx->idB_len,
            &ctx->RB)) {
        fprintf(stderr, "setup failed: initiator precompute S1\n");
        return 0;
    }

    {
        sm2_affine_t RBc = ctx->RB;
        if (!sm2_kex_responder_compute_key_mont(
                &RBc, ctx->KB, 32, ctx->S2, ctx->SB, ctx->S1,
                ctx->dB, &ctx->PB, ctx->idB, ctx->idB_len,
                ctx->rB, &ctx->RA, &ctx->PA, ctx->idA, ctx->idA_len)) {
            fprintf(stderr, "setup failed: responder precompute S2\n");
            return 0;
        }
    }

    sink_mem(&ctx->PA, sizeof(ctx->PA));
    sink_mem(&ctx->PB, sizeof(ctx->PB));
    sink_mem(&ctx->RA, sizeof(ctx->RA));
    sink_mem(&ctx->RB, sizeof(ctx->RB));
    sink_mem(ctx->S1, sizeof(ctx->S1));
    sink_mem(ctx->S2, sizeof(ctx->S2));

    return 1;
}

static void warmup_genra(const perf_ctx_t *ctx, unsigned long n) {
    sm2_affine_t R;
    for (unsigned long i = 0; i < n; i++) {
        sm2_kex_initiator_gen_RA_mont(&R, ctx->rA);
        sink_mem(&R, sizeof(R));
    }
}

static void run_genra(const perf_ctx_t *ctx, unsigned long iters) {
    sm2_affine_t R;
    for (unsigned long i = 0; i < iters; i++) {
        if (!sm2_kex_initiator_gen_RA_mont(&R, ctx->rA)) {
            fprintf(stderr, "genra failed at iter %lu\n", i);
            exit(1);
        }
        sink_mem(&R, sizeof(R));
    }
}

static void warmup_scalar(const perf_ctx_t *ctx, unsigned long n) {
    sm2_affine_t Gm;
    sm2_jacobian_t RmJ;
    sm2_affine_to_mont(&Gm, &ctx->G);
    for (unsigned long i = 0; i < n; i++) {
        sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&RmJ, &Gm, ctx->dA);
        sink_mem(&RmJ, sizeof(RmJ));
    }
}

static void run_scalar(const perf_ctx_t *ctx, unsigned long iters) {
    sm2_affine_t Gm;
    sm2_jacobian_t RmJ;
    sm2_affine_to_mont(&Gm, &ctx->G);
    for (unsigned long i = 0; i < iters; i++) {
        sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&RmJ, &Gm, ctx->dA);
        sink_mem(&RmJ, sizeof(RmJ));
    }
}

static void warmup_responder(const perf_ctx_t *ctx, unsigned long n) {
    sm2_affine_t RBc;
    uint8_t KB[32], S2[32], SB[32];
    for (unsigned long i = 0; i < n; i++) {
        RBc = ctx->RB;
        sm2_kex_responder_compute_key_mont(
            &RBc, KB, 32, S2, SB, ctx->S1,
            ctx->dB, &ctx->PB, ctx->idB, ctx->idB_len,
            ctx->rB, &ctx->RA, &ctx->PA, ctx->idA, ctx->idA_len);
        sink_mem(&RBc, sizeof(RBc));
        sink_mem(KB, sizeof(KB));
        sink_mem(S2, sizeof(S2));
        sink_mem(SB, sizeof(SB));
    }
}

static void run_responder(const perf_ctx_t *ctx, unsigned long iters) {
    sm2_affine_t RBc;
    uint8_t KB[32], S2[32], SB[32];
    for (unsigned long i = 0; i < iters; i++) {
        RBc = ctx->RB;
        if (!sm2_kex_responder_compute_key_mont(
                &RBc, KB, 32, S2, SB, ctx->S1,
                ctx->dB, &ctx->PB, ctx->idB, ctx->idB_len,
                ctx->rB, &ctx->RA, &ctx->PA, ctx->idA, ctx->idA_len)) {
            fprintf(stderr, "responder failed at iter %lu\n", i);
            exit(1);
        }
        sink_mem(&RBc, sizeof(RBc));
        sink_mem(KB, sizeof(KB));
        sink_mem(S2, sizeof(S2));
        sink_mem(SB, sizeof(SB));
    }
}

static void warmup_initiator(const perf_ctx_t *ctx, unsigned long n) {
    uint8_t KA[32], S1[32], SA[32];
    for (unsigned long i = 0; i < n; i++) {
        sm2_kex_initiator_compute_key_mont(
            KA, 32, S1, SA, ctx->S2,
            ctx->dA, &ctx->PA, ctx->idA, ctx->idA_len,
            ctx->rA, &ctx->RA,
            &ctx->PB, ctx->idB, ctx->idB_len,
            &ctx->RB);
        sink_mem(KA, sizeof(KA));
        sink_mem(S1, sizeof(S1));
        sink_mem(SA, sizeof(SA));
    }
}

static void run_initiator(const perf_ctx_t *ctx, unsigned long iters) {
    uint8_t KA[32], S1[32], SA[32];
    for (unsigned long i = 0; i < iters; i++) {
        if (!sm2_kex_initiator_compute_key_mont(
                KA, 32, S1, SA, ctx->S2,
                ctx->dA, &ctx->PA, ctx->idA, ctx->idA_len,
                ctx->rA, &ctx->RA,
                &ctx->PB, ctx->idB, ctx->idB_len,
                &ctx->RB)) {
            fprintf(stderr, "initiator failed at iter %lu\n", i);
            exit(1);
        }
        sink_mem(KA, sizeof(KA));
        sink_mem(S1, sizeof(S1));
        sink_mem(SA, sizeof(SA));
    }
}

static void warmup_full(const perf_ctx_t *ctx, unsigned long n) {
    sm2_affine_t RBc;
    uint8_t KB[32], S2[32], SB[32];
    uint8_t KA[32], S1[32], SA[32];

    for (unsigned long i = 0; i < n; i++) {
        RBc = ctx->RB;
        sm2_kex_responder_compute_key_mont(
            &RBc, KB, 32, S2, SB, ctx->S1,
            ctx->dB, &ctx->PB, ctx->idB, ctx->idB_len,
            ctx->rB, &ctx->RA, &ctx->PA, ctx->idA, ctx->idA_len);

        sm2_kex_initiator_compute_key_mont(
            KA, 32, S1, SA, S2,
            ctx->dA, &ctx->PA, ctx->idA, ctx->idA_len,
            ctx->rA, &ctx->RA,
            &ctx->PB, ctx->idB, ctx->idB_len,
            &ctx->RB);

        sink_mem(&RBc, sizeof(RBc));
        sink_mem(KB, sizeof(KB));
        sink_mem(S2, sizeof(S2));
        sink_mem(SB, sizeof(SB));
        sink_mem(KA, sizeof(KA));
        sink_mem(S1, sizeof(S1));
        sink_mem(SA, sizeof(SA));
    }
}

static void run_full(const perf_ctx_t *ctx, unsigned long iters) {
    sm2_affine_t RBc;
    uint8_t KB[32], S2[32], SB[32];
    uint8_t KA[32], S1[32], SA[32];

    for (unsigned long i = 0; i < iters; i++) {
        RBc = ctx->RB;
        if (!sm2_kex_responder_compute_key_mont(
                &RBc, KB, 32, S2, SB, ctx->S1,
                ctx->dB, &ctx->PB, ctx->idB, ctx->idB_len,
                ctx->rB, &ctx->RA, &ctx->PA, ctx->idA, ctx->idA_len)) {
            fprintf(stderr, "full/responder failed at iter %lu\n", i);
            exit(1);
        }

        if (!sm2_kex_initiator_compute_key_mont(
                KA, 32, S1, SA, S2,
                ctx->dA, &ctx->PA, ctx->idA, ctx->idA_len,
                ctx->rA, &ctx->RA,
                &ctx->PB, ctx->idB, ctx->idB_len,
                &ctx->RB)) {
            fprintf(stderr, "full/initiator failed at iter %lu\n", i);
            exit(1);
        }

        sink_mem(&RBc, sizeof(RBc));
        sink_mem(KB, sizeof(KB));
        sink_mem(S2, sizeof(S2));
        sink_mem(SB, sizeof(SB));
        sink_mem(KA, sizeof(KA));
        sink_mem(S1, sizeof(S1));
        sink_mem(SA, sizeof(SA));
    }
}

static int run_mode(const perf_ctx_t *ctx, const char *mode, unsigned long iters) {
    const unsigned long warmup = 16UL;

    if (strcmp(mode, "genra") == 0) {
        printf("[mode] genra, iters=%lu\n", iters);
        warmup_genra(ctx, warmup);
        run_genra(ctx, iters);
        return 1;
    }

    if (strcmp(mode, "scalar") == 0) {
        printf("[mode] scalar, iters=%lu\n", iters);
        warmup_scalar(ctx, warmup);
        run_scalar(ctx, iters);
        return 1;
    }

    if (strcmp(mode, "responder") == 0) {
        printf("[mode] responder, iters=%lu\n", iters);
        warmup_responder(ctx, warmup);
        run_responder(ctx, iters);
        return 1;
    }

    if (strcmp(mode, "initiator") == 0) {
        printf("[mode] initiator, iters=%lu\n", iters);
        warmup_initiator(ctx, warmup);
        run_initiator(ctx, iters);
        return 1;
    }

    if (strcmp(mode, "full") == 0) {
        printf("[mode] full, iters=%lu\n", iters);
        warmup_full(ctx, warmup);
        run_full(ctx, iters);
        return 1;
    }

    if (strcmp(mode, "all") == 0) {
        printf("[mode] all, per-mode iters=%lu\n", iters);

        warmup_genra(ctx, warmup);
        run_genra(ctx, iters);

        warmup_scalar(ctx, warmup);
        run_scalar(ctx, iters);

        warmup_responder(ctx, warmup);
        run_responder(ctx, iters);

        warmup_initiator(ctx, warmup);
        run_initiator(ctx, iters);

        warmup_full(ctx, warmup);
        run_full(ctx, iters);

        return 1;
    }

    return 0;
}

int main(int argc, char **argv) {
    const char *mode = "all";
    unsigned long iters = DEFAULT_PERF_ITERS;
    perf_ctx_t ctx;
    char *endp = NULL;

    if (argc >= 2) {
        if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
            usage(argv[0]);
            return 0;
        }
        mode = argv[1];
    }

    if (argc >= 3) {
        iters = strtoul(argv[2], &endp, 10);
        if (endp == argv[2] || *endp != '\0' || iters == 0) {
            fprintf(stderr, "invalid iteration count: %s\n", argv[2]);
            return 1;
        }
    }

    if (!setup_perf_ctx(&ctx)) {
        return 1;
    }

    if (!run_mode(&ctx, mode, iters)) {
        fprintf(stderr, "unknown mode: %s\n\n", mode);
        usage(argv[0]);
        return 1;
    }

    printf("[sink] %u\n", bench_sink);
    return 0;
}
