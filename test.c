/*
 * benchmark.c
 * ---------------------------------------
 * SM2 constant-time performance benchmark
 * SchemeB window scalar multiplication
 * RISC-V rdcycle based measurement
 *
 * This benchmark is:
 *  - reproducible
 *  - multi-round averaged
 *  - optimizer-safe
 *
 * Author: (you)
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#include "fn.h"
#include "fp.h"
#include "sm2_curve.h"
#include "sm2_scalar.h"
#include "sm2_kex.h"

/* ============================================================
 * RISC-V cycle counter
 * ============================================================ */
static inline uint64_t rdcycle(void) {
    uint64_t c;
    asm volatile ("rdcycle %0" : "=r"(c));
    return c;
}

/* ============================================================
 * Anti-optimization sink
 * ============================================================ */
volatile uint8_t bench_sink;

static inline void sink_mem(const void *p, size_t n) {
    const volatile uint8_t *b = (const volatile uint8_t *)p;
    for (size_t i = 0; i < n; i++) {
        bench_sink ^= b[i];
    }
}

/* ============================================================
 * Benchmark parameters
 * ============================================================ */
#define REPEAT_FP     5
#define REPEAT_MUL    5
#define REPEAT_KEX    5

#define ITER_FP       20000
#define ITER_MUL      200
#define ITER_KEX      50

/* ============================================================
 * Statistics helper
 * ============================================================ */
static void print_stats(const char *name, uint64_t *v, int n) {
    uint64_t sum = 0, min = v[0], max = v[0];
    for (int i = 0; i < n; i++) {
        sum += v[i];
        if (v[i] < min) min = v[i];
        if (v[i] > max) max = v[i];
    }
    printf("%-30s : mean=%" PRIu64
           ", min=%" PRIu64
           ", max=%" PRIu64 "\n",
           name, sum / n, min, max);
}

/* ============================================================
 * Deterministic test data
 * ============================================================ */
static void fixed_scalar(uint8_t k[32], uint8_t v) {
    memset(k, v, 32);
    k[31] |= 1; /* ensure non-zero & odd */
}

/* ============================================================
 * Main benchmark
 * ============================================================ */
int main(void) {
    printf("SM2 Benchmark (strict, reproducible)\n");
    printf("===================================\n");
    printf("iters: fp=%d, mul=%d, kex=%d\n\n",
           ITER_FP, ITER_MUL, ITER_KEX);

    /* ---------- common data ---------- */
    fp_t a, b, r;
    sm2_affine_t G, PA, PB, RA, RB;
    sm2_jacobian_t J;

    uint8_t dA[32], dB[32], rA[32], rB[32];
    uint8_t KA[32], KB[32];
    uint8_t S1[32], S2[32], SA[32], SB[32];

    const uint8_t idA[] = "Alice";
    const uint8_t idB[] = "Bob";

    sm2_get_base_affine(&G);

    fixed_scalar(dA, 0x11);
    fixed_scalar(dB, 0x22);
    fixed_scalar(rA, 0x33);
    fixed_scalar(rB, 0x44);
    
    
    
    /* ========================================================
 * fp_reduce benchmark (512->256 reduction, chained)
 *   assume: void fp_reduce(fp_t *r, const uint64_t T[8]);
 * ======================================================== */
uint64_t red_cycles[REPEAT_FP];

uint64_t T[8];
for (int i = 0; i < 8; i++) {
    T[i] = 0x1111111111111111ULL * (uint64_t)(i + 1);
}

fp_t rr;
memset(&rr, 0, sizeof(rr));

for (int r0 = 0; r0 < REPEAT_FP; r0++) {
    /* warm-up */
    for (int i = 0; i < 200; i++) {
        fp_reduce(&rr, T);
        T[0] ^= rr.v[0];
    }

    uint64_t start = rdcycle();
    for (int i = 0; i < ITER_FP; i++) {
        fp_reduce(&rr, T);

        /* chained dependency: next input depends on output */
        T[0] += rr.v[0];
        T[1] ^= rr.v[1];
        T[2] += rr.v[2];
        T[3] ^= rr.v[3];
    }
    uint64_t end = rdcycle();

    red_cycles[r0] = (end - start) / (uint64_t)ITER_FP;
    sink_mem(&rr, sizeof(rr));
    sink_mem(T, sizeof(T));
}

print_stats("fp_reduce (chained)", red_cycles, REPEAT_FP);

    

    /* ========================================================
     * fp_mul benchmark
     * ======================================================== */
/* ========================================================
     * fp_mul benchmark (Strict & High Precision)
     * ======================================================== */
    uint64_t fp_cycles[REPEAT_FP];
    
    // 准备链式数据，防止编译器优化
    fp_t t_data;
    memset(&t_data, 0x12, sizeof(fp_t)); 

    for (int r0 = 0; r0 < REPEAT_FP; r0++) {
        // 1. 测量空循环开销 (Calibration)
        uint64_t c_start = rdcycle();
        for (int i = 0; i < ITER_FP; i++) {
            // 空操作，但保留 sink 防止循环被完全优化
            sink_mem(&t_data, sizeof(t_data));
        }
        uint64_t c_end = rdcycle();
        uint64_t empty_loop_total = c_end - c_start;

        // 2. 测量实际 fp_mul (带链式依赖)
        uint64_t start = rdcycle();
        for (int i = 0; i < ITER_FP; i++) {
            // 让结果参与下一次运算，强制顺序执行
            fp_mul(&t_data, &t_data, &t_data); 
            // 每次计算后手动 sink 一次
            if (i % 100 == 0) sink_mem(&t_data, sizeof(t_data));
        }
        uint64_t end = rdcycle();

        // 计算纯净的函数周期
        uint64_t total_cycles = (end - start);
        if (total_cycles > empty_loop_total) {
            fp_cycles[r0] = (total_cycles - empty_loop_total) / ITER_FP;
        } else {
            fp_cycles[r0] = total_cycles / ITER_FP; // 防止意外
        }
    }

    print_stats("fp_mul (Precise)", fp_cycles, REPEAT_FP);
    
    /* ========================================================
 * fp_sqr benchmark (dependency-chained)
 * ======================================================== */
uint64_t sqr_cycles[REPEAT_FP];

fp_t xs;
memset(&xs, 0x56, sizeof(xs));
xs.v[0] |= 1;

for (int r0 = 0; r0 < REPEAT_FP; r0++) {
    /* warm-up */
    for (int i = 0; i < 200; i++) {
        fp_sqr(&xs, &xs);
    }

    uint64_t start = rdcycle();
    for (int i = 0; i < ITER_FP; i++) {
        fp_sqr(&xs, &xs);  /* chained */
    }
    uint64_t end = rdcycle();

    sqr_cycles[r0] = (end - start) / (uint64_t)ITER_FP;
    sink_mem(&xs, sizeof(xs));
}

print_stats("fp_sqr (chained)", sqr_cycles, REPEAT_FP);


   /* ========================================================
 * fp_inv benchmark (Fermat / whatever fp_inv uses)
 * ======================================================== */
#define ITER_INV  200   /* 先用 200；如果太快可加到 500/1000 */
uint64_t inv_cycles[REPEAT_FP];

fp_t xi, ri;
memset(&xi, 0x9a, sizeof(xi));
xi.v[0] |= 1; /* non-zero */

for (int r0 = 0; r0 < REPEAT_FP; r0++) {
    /* warm-up */
    for (int i = 0; i < 10; i++) {
        fp_inv(&ri, &xi);
        sink_mem(&ri, sizeof(ri));
    }

    uint64_t start = rdcycle();
    for (int i = 0; i < ITER_INV; i++) {
        fp_inv(&ri, &xi);

        /* 轻微扰动输入，避免编译器/模拟器走完全相同路径 */
        fp_add(&xi, &xi, &ri);
        xi.v[0] |= 1;
    }
    uint64_t end = rdcycle();

    inv_cycles[r0] = (end - start) / (uint64_t)ITER_INV;
    sink_mem(&ri, sizeof(ri));
}

print_stats("fp_inv (avg)", inv_cycles, REPEAT_FP);

    
    
    /* ========================================================
     * Scalar multiplication (SchemeB)
     * ======================================================== */
    uint64_t mul_cycles[REPEAT_MUL];

    for (int r0 = 0; r0 < REPEAT_MUL; r0++) {
        for (int i = 0; i < 10; i++) {
            sm2_scalar_mul_window_ct_schemeB(&J, &G, dA);
        }

        uint64_t start = rdcycle();
        for (int i = 0; i < ITER_MUL; i++) {
            sm2_scalar_mul_window_ct_schemeB(&J, &G, dA);
            sink_mem(&J, sizeof(J));
        }
        uint64_t end = rdcycle();

        mul_cycles[r0] = (end - start) / ITER_MUL;
    }

    print_stats("ScalarMul SchemeB ([d]G)", mul_cycles, REPEAT_MUL);

    /* ========================================================
     * KEX: Initiator generates RA
     * ======================================================== */
    uint64_t kex_gen_cycles[REPEAT_KEX];

    for (int r0 = 0; r0 < REPEAT_KEX; r0++) {
        for (int i = 0; i < 5; i++) {
            sm2_kex_initiator_gen_RA(&RA, rA);
        }

        uint64_t start = rdcycle();
        for (int i = 0; i < ITER_KEX; i++) {
            sm2_kex_initiator_gen_RA(&RA, rA);
            sink_mem(&RA, sizeof(RA));
        }
        uint64_t end = rdcycle();

        kex_gen_cycles[r0] = (end - start) / ITER_KEX;
    }

    print_stats("KEX Initiator Gen RA", kex_gen_cycles, REPEAT_KEX);

    /* ========================================================
     * Prepare long-term public keys
     * ======================================================== */
    sm2_scalar_mul_window_ct_schemeB(&J, &G, dA);
    sm2_jacobian_to_affine(&PA, &J);

    sm2_scalar_mul_window_ct_schemeB(&J, &G, dB);
    sm2_jacobian_to_affine(&PB, &J);

    sm2_kex_initiator_gen_RA(&RA, rA);

    /* ========================================================
     * KEX: Responder full path
     * ======================================================== */
    uint64_t kex_resp_cycles[REPEAT_KEX];

    for (int r0 = 0; r0 < REPEAT_KEX; r0++) {
        for (int i = 0; i < 3; i++) {
            sm2_kex_responder_compute_key(
                &RB, KB, 32, S2, SB, S1,
                dB, &PB, idB, sizeof(idB) - 1,
                rB, &RA, &PA, idA, sizeof(idA) - 1
            );
        }

        uint64_t start = rdcycle();
        for (int i = 0; i < ITER_KEX; i++) {
            sm2_kex_responder_compute_key(
                &RB, KB, 32, S2, SB, S1,
                dB, &PB, idB, sizeof(idB) - 1,
                rB, &RA, &PA, idA, sizeof(idA) - 1
            );
            sink_mem(KB, 32);
        }
        uint64_t end = rdcycle();

        kex_resp_cycles[r0] = (end - start) / ITER_KEX;
    }

    print_stats("KEX Responder (full path)", kex_resp_cycles, REPEAT_KEX);

    /* ========================================================
     * KEX: Initiator verify path
     * ======================================================== */
    uint64_t kex_init_cycles[REPEAT_KEX];

    for (int r0 = 0; r0 < REPEAT_KEX; r0++) {
        for (int i = 0; i < 3; i++) {
            sm2_kex_initiator_compute_key(
                KA, 32, S1, SA, S2,
                dA, &PA, idA, sizeof(idA) - 1,
                rA, &RA,
                &PB, idB, sizeof(idB) - 1,
                &RB
            );
        }

        uint64_t start = rdcycle();
        for (int i = 0; i < ITER_KEX; i++) {
            sm2_kex_initiator_compute_key(
                KA, 32, S1, SA, S2,
                dA, &PA, idA, sizeof(idA) - 1,
                rA, &RA,
                &PB, idB, sizeof(idB) - 1,
                &RB
            );
            sink_mem(KA, 32);
        }
        uint64_t end = rdcycle();

        kex_init_cycles[r0] = (end - start) / ITER_KEX;
    }

    print_stats("KEX Initiator (full verify)", kex_init_cycles, REPEAT_KEX);

    printf("\n[sink] %u\n", bench_sink);
    return 0;
}
