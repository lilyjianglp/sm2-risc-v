/*
 * test.c  (原 benchmark.c)
 * ---------------------------------------
 * SM2 constant-time performance benchmark
 * SchemeB window scalar multiplication
 * RISC-V rdcycle based measurement
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


static int bytes_equal(const uint8_t* a, const uint8_t* b, size_t len) {
    return memcmp(a, b, len) == 0;
}

static void print_hex(const char* label, const uint8_t* data, size_t len) {
    printf("%s = ", label);
    for (size_t i = 0; i < len; i++) printf("%02x", data[i]);
    printf("\n");
}

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
    for (size_t i = 0; i < n; i++) bench_sink ^= b[i];
}

/* ============================================================
 * Deterministic PRNG (xorshift64)
 * ============================================================ */
static uint64_t g_state = 0xC0FFEE123456789AULL;

static inline uint64_t xorshift64(void) {
    uint64_t x = g_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    g_state = x;
    return x;
}

static inline void rand_bytes(uint8_t *out, size_t n) {
    for (size_t i = 0; i < n; ) {
        uint64_t r = xorshift64();
        for (int j = 0; j < 8 && i < n; j++, i++)
            out[i] = (uint8_t)(r >> (56 - 8 * j));
    }
}

static inline int bytes_is_zero32(const uint8_t a[32]) {
    uint8_t v = 0;
    for (int i = 0; i < 32; i++) v |= a[i];
    return v == 0;
}

static inline void rand_scalar_mod_n(uint8_t out[32]) {
    fn_t t;
    uint8_t raw[32];
    do {
        rand_bytes(raw, 32);
        fn_from_be(&t, raw);
        fn_to_be(out, &t);
        sink_mem(out, 32);
    } while (bytes_is_zero32(out));
}

/* ============================================================
 * Benchmark parameters
 * ============================================================ */
#ifndef REPEAT_FP
#define REPEAT_FP   5
#endif

#ifndef REPEAT_MUL
#define REPEAT_MUL  5
#endif

#ifndef REPEAT_KEX
#define REPEAT_KEX  5
#endif

#ifndef REPEAT_FN
#define REPEAT_FN   5
#endif

#ifndef REPEAT_INV
#define REPEAT_INV  9
#endif

#ifndef ITER_FP
#define ITER_FP     100000
#endif

#ifndef ITER_MUL
#define ITER_MUL    1000
#endif

#ifndef ITER_KEX
#define ITER_KEX    200
#endif

#ifndef ITER_INV
#define ITER_INV    10000
#endif

#ifndef ITER_FN
#define ITER_FN     20000
#endif

#ifdef USE_RV64_ASM_MUL
#define FP_MONT_MUL_LABEL "fp_mont_mul (ASM)"
#else
#define FP_MONT_MUL_LABEL "fp_mont_mul (C)"
#endif

#ifdef USE_RV64_ASM_SQR
#define FP_MONT_SQR_LABEL "fp_mont_sqr (ASM)"
#else
#define FP_MONT_SQR_LABEL "fp_mont_sqr (C)"
#endif

#ifdef USE_PORNIN_FULL_INV_ASM
#define FP_MONT_INV_LABEL "fp_mont_inv (Pornin full ASM)"
#elif defined(USE_PORNIN_INV) && defined(USE_PORNIN_INNER31_ASM) && defined(USE_PORNIN_UPDATE_UV_ASM)
#define FP_MONT_INV_LABEL "fp_mont_inv (Pornin C+ASM)"
#elif defined(USE_PORNIN_INV) && defined(USE_PORNIN_INNER31_ASM)
#define FP_MONT_INV_LABEL "fp_mont_inv (Pornin inner31 ASM)"
#elif defined(USE_PORNIN_INV)
#define FP_MONT_INV_LABEL "fp_mont_inv (Pornin C)"
#else
#define FP_MONT_INV_LABEL "fp_mont_inv (C)"
#endif

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
    printf("%-38s : mean=%" PRIu64 ", min=%" PRIu64 ", max=%" PRIu64 "\n",
           name, sum / (uint64_t)n, min, max);
}

static void fixed_scalar(uint8_t k[32], uint8_t v) {
    memset(k, v, 32);
    k[31] |= 1;
}

/* ============================================================
 * Correctness helpers
 * ============================================================ */
static int fp_equal_local(const fp_t *a, const fp_t *b) { return fp_is_equal(a, b); }

static int point_equal(const sm2_affine_t *P, const sm2_affine_t *Q) {
    if (P->infinity != Q->infinity) return 0;
    if (P->infinity) return 1;
    return fp_equal_local(&P->x, &Q->x) && fp_equal_local(&P->y, &Q->y);
}

static void point_to_bytes(uint8_t out[64], const sm2_affine_t *P) {
    fp_to_bytes(out, &P->x);
    fp_to_bytes(out + 32, &P->y);
}

static void print_point(const char *label, const sm2_affine_t *P) {
    uint8_t buf[64];
    if (P->infinity) { printf("%s = INF\n", label); return; }
    point_to_bytes(buf, P);
    printf("%s.x = ", label);
    for (int i = 0; i < 32; i++) printf("%02x", buf[i]);
    printf("\n");
    printf("%s.y = ", label);
    for (int i = 32; i < 64; i++) printf("%02x", buf[i]);
    printf("\n");
}

static void scalar_set_zero(uint8_t k[32])           { memset(k, 0, 32); }

static void scalar_set_u64(uint8_t k[32], uint64_t x) {
    memset(k, 0, 32);
    for (int i = 0; i < 8; i++) { k[31 - i] = (uint8_t)(x & 0xFF); x >>= 8; }
}

static void scalar_set_pow2(uint8_t k[32], int bit) {
    memset(k, 0, 32);
    k[31 - (bit >> 3)] = (uint8_t)(1U << (bit & 7));
}

static int is_zero_scalar(const uint8_t k[32]) {
    uint8_t x = 0;
    for (int i = 0; i < 32; i++) x |= k[i];
    return x == 0;
}

static uint64_t g_test_rng_state = 0x123456789abcdef0ULL;

static uint64_t test_xorshift64(void) {
    uint64_t x = g_test_rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    g_test_rng_state = x;
    return x;
}

static void scalar_random_nonzero_test(uint8_t k[32]) {
    do {
        for (int i = 0; i < 4; i++) {
            uint64_t w = test_xorshift64();
            for (int j = 0; j < 8; j++)
                k[i * 8 + j] = (uint8_t)(w >> (56 - 8 * j));
        }
        k[0] &= 0x7F;
    } while (is_zero_scalar(k));
}

static void scalar_mul_normal_affine(sm2_affine_t *R, const sm2_affine_t *P, const uint8_t k[32]) {
    sm2_jacobian_t J;
    sm2_scalar_mul_window_ct_schemeB(&J, P, k);
    sm2_jacobian_to_affine(R, &J);
}

#ifndef USE_FP_MONT
static void scalar_mul_normal_ladder_affine(sm2_affine_t *R, const sm2_affine_t *P, const uint8_t k[32]) {
    sm2_jacobian_t J;
    sm2_scalar_mul_ladder_ct(&J, P, k);
    sm2_jacobian_to_affine(R, &J);
}
#endif /* !USE_FP_MONT */

#ifdef USE_FP_MONT
static void scalar_mul_mont_affine(sm2_affine_t *R, const sm2_affine_t *P, const uint8_t k[32]) {
    sm2_scalar_mul_window_ct_schemeB_mont_to_affine(R, P, k);
}

static void derive_public_key_mont(sm2_affine_t *Ppub, const uint8_t d[32]) {
    sm2_affine_t G; sm2_get_base_affine(&G);
    scalar_mul_mont_affine(Ppub, &G, d);
}
#endif /* USE_FP_MONT */

static void derive_public_key_normal(sm2_affine_t *Ppub, const uint8_t d[32]) {
    sm2_affine_t G; sm2_get_base_affine(&G);
    scalar_mul_normal_affine(Ppub, &G, d);
}

/* ── KEX 正确性辅助 ── */
static int run_one_kex_normal(uint8_t K[32], uint8_t S1[32], uint8_t S2[32],
                              const uint8_t dA[32], const uint8_t dB[32],
                              const uint8_t rA[32], const uint8_t rB[32],
                              const sm2_affine_t *PA, const sm2_affine_t *PB,
                              const uint8_t *idA, size_t idA_len,
                              const uint8_t *idB, size_t idB_len)
{
    sm2_affine_t RA, RB;
    uint8_t KA1[32], KA2[32], KB[32], SA1[32], SA2[32], SB[32];
    uint8_t S1l[32], S2l[32];

    if (!sm2_kex_initiator_gen_RA(&RA, rA)) return 0;
    if (!sm2_kex_initiator_gen_RA(&RB, rB)) return 0;

    if (!sm2_kex_initiator_compute_key(KA1,32,S1l,SA1,NULL,
            dA,PA,idA,idA_len,rA,&RA,PB,idB,idB_len,&RB)) return 0;

    { sm2_affine_t RBc = RB;
      if (!sm2_kex_responder_compute_key(&RBc,KB,32,S2l,SB,S1l,
              dB,PB,idB,idB_len,rB,&RA,PA,idA,idA_len)) return 0; }

    if (!sm2_kex_initiator_compute_key(KA2,32,S1,SA2,S2l,
            dA,PA,idA,idA_len,rA,&RA,PB,idB,idB_len,&RB)) return 0;

    if (memcmp(KA1,KB,32)||memcmp(KA2,KB,32)||
        memcmp(SA1,S2l,32)||memcmp(S1,S1l,32)) return 0;
    memcpy(K,KB,32); memcpy(S2,S2l,32);
    return 1;
}

#ifdef USE_FP_MONT
static int run_one_kex_mont(uint8_t K[32], uint8_t S1[32], uint8_t S2[32],
                            const uint8_t dA[32], const uint8_t dB[32],
                            const uint8_t rA[32], const uint8_t rB[32],
                            const sm2_affine_t *PA, const sm2_affine_t *PB,
                            const uint8_t *idA, size_t idA_len,
                            const uint8_t *idB, size_t idB_len)
{
    sm2_affine_t RA, RB;
    uint8_t KA1[32], KA2[32], KB[32], SA1[32], SA2[32], SB[32];
    uint8_t S1l[32], S2l[32];

    if (!sm2_kex_initiator_gen_RA_mont(&RA, rA)) return 0;
    if (!sm2_kex_initiator_gen_RA_mont(&RB, rB)) return 0;

    if (!sm2_kex_initiator_compute_key_mont(KA1,32,S1l,SA1,NULL,
            dA,PA,idA,idA_len,rA,&RA,PB,idB,idB_len,&RB)) return 0;

    { sm2_affine_t RBc = RB;
      if (!sm2_kex_responder_compute_key_mont(&RBc,KB,32,S2l,SB,S1l,
              dB,PB,idB,idB_len,rB,&RA,PA,idA,idA_len)) return 0; }

    if (!sm2_kex_initiator_compute_key_mont(KA2,32,S1,SA2,S2l,
            dA,PA,idA,idA_len,rA,&RA,PB,idB,idB_len,&RB)) return 0;

    if (memcmp(KA1,KB,32)||memcmp(KA2,KB,32)||
        memcmp(SA1,S2l,32)||memcmp(S1,S1l,32)) return 0;
    memcpy(K,KB,32); memcpy(S2,S2l,32);
    return 1;
}

/* ── 域转换正确性 ── */
static int test_domain_conversion(const sm2_affine_t *P, const char *name) {
    sm2_affine_t Pm, Pback, A_from_J, A_mont, A_mont_back;
    sm2_jacobian_t Jm;

    sm2_affine_to_mont(&Pm, P);
    sm2_affine_from_mont(&Pback, &Pm);
    if (!point_equal(P, &Pback)) {
        printf("FAIL: domain conversion affine round-trip (%s)\n", name); return 0;
    }

    sm2_affine_to_jacobian_mont(&Jm, P);
    sm2_jacobian_to_affine_mont(&A_from_J, &Jm);
    if (!point_equal(P, &A_from_J)) {
        printf("FAIL: domain conversion jacobian->affine round-trip (%s)\n", name); return 0;
    }

    sm2_affine_to_mont(&Pm, P);
    sm2_affine_mont_to_jacobian_mont(&Jm, &Pm);
    sm2_jacobian_to_affine_mont(&A_mont, &Jm);
    sm2_affine_from_mont(&A_mont_back, &A_mont);
    if (!point_equal(P, &A_mont_back)) {
        printf("FAIL: domain conversion full Montgomery round-trip (%s)\n", name); return 0;
    }

    printf("PASS: domain conversion (%s)\n", name);
    return 1;
}

static int test_curve_ops_diff(const sm2_affine_t *P, const sm2_affine_t *Q, const char *name) {
    sm2_jacobian_t Jp, Jq, Jp_m, Jq_m;
    sm2_affine_t Rn, Rm;

    sm2_affine_to_jacobian(&Jp, P); sm2_affine_to_jacobian(&Jq, Q);
    sm2_affine_to_jacobian_mont(&Jp_m, P); sm2_affine_to_jacobian_mont(&Jq_m, Q);

    { sm2_jacobian_t Tn, Tm;
      sm2_double_jm_a_minus3(&Tn, &Jp); sm2_double_jm_a_minus3_mont(&Tm, &Jp_m);
      sm2_jacobian_to_affine(&Rn, &Tn); sm2_jacobian_to_affine_mont(&Rm, &Tm);
      if (!point_equal(&Rn, &Rm)) { printf("FAIL: curve op double (%s)\n", name); return 0; } }

    { sm2_jacobian_t Tn, Tm; sm2_affine_t Qm;
      sm2_affine_to_mont(&Qm, Q);
      sm2_add_ja(&Tn, &Jp, Q); sm2_add_ja_mont(&Tm, &Jp_m, &Qm);
      sm2_jacobian_to_affine(&Rn, &Tn); sm2_jacobian_to_affine_mont(&Rm, &Tm);
      if (!point_equal(&Rn, &Rm)) { printf("FAIL: curve op add_ja (%s)\n", name); return 0; } }

    { sm2_jacobian_t Tn, Tm;
      sm2_add_jj(&Tn, &Jp, &Jq); sm2_add_jj_mont(&Tm, &Jp_m, &Jq_m);
      sm2_jacobian_to_affine(&Rn, &Tn); sm2_jacobian_to_affine_mont(&Rm, &Tm);
      if (!point_equal(&Rn, &Rm)) { printf("FAIL: curve op add_jj (%s)\n", name); return 0; } }

    printf("PASS: curve ops diff (%s)\n", name);
    return 1;
}

static int run_scalar_diff_case(const char *name, const sm2_affine_t *P,
                                const uint8_t k[32], int verbose) {
    sm2_affine_t Rn, Rm;
    scalar_mul_normal_affine(&Rn, P, k);
    scalar_mul_mont_affine(&Rm, P, k);
    if (!point_equal(&Rn, &Rm)) {
        printf("FAIL: scalar diff (%s)\n", name);
        if (verbose) { print_hex("k",k,32); print_point("normal",&Rn); print_point("mont  ",&Rm); }
        return 0;
    }
    printf("PASS: scalar diff (%s)\n", name);
    return 1;
}

#ifdef USE_FP_MONT
static int run_fixedbase_special_case(const char *name, const uint8_t k[32]) {
    sm2_affine_t G, Rref, Rmont;
    sm2_get_base_affine(&G);
    scalar_mul_normal_affine(&Rref, &G, k);
    scalar_mul_mont_affine(&Rmont, &G, k);
    if (!point_equal(&Rref, &Rmont)) {
        printf("FAIL: fixed-base special (%s)\n", name);
        print_hex("k",k,32); print_point("ladder_ref",&Rref); print_point("mont_fixed",&Rmont);
        return 0;
    }
    printf("PASS: fixed-base special (%s)\n", name);
    return 1;
}
#endif /* USE_FP_MONT */

static int run_comprehensive_correctness_suite(void) {
    static const uint8_t tdA[32] = {
        0x12,0x34,0x56,0x78,0x9a,0xbc,0xde,0xf0,0x11,0x22,0x33,0x44,0x55,0x66,0x77,0x88,
        0x99,0xaa,0xbb,0xcc,0xdd,0xee,0xf1,0x02,0x13,0x24,0x35,0x46,0x57,0x68,0x79,0x8a };
    static const uint8_t tdB[32] = {
        0x23,0x45,0x67,0x89,0xab,0xcd,0xef,0x01,0x10,0x32,0x54,0x76,0x98,0xba,0xdc,0xfe,
        0x12,0x34,0x56,0x78,0x87,0x65,0x43,0x21,0x0f,0x1e,0x2d,0x3c,0x4b,0x5a,0x69,0x78 };
    static const uint8_t trA[32] = {
        0x34,0x56,0x78,0x9a,0xbc,0xde,0xf0,0x12,0x21,0x43,0x65,0x87,0xa9,0xcb,0xed,0x0f,
        0x13,0x57,0x9b,0xdf,0x24,0x68,0xac,0xe0,0x10,0x20,0x30,0x40,0x50,0x60,0x70,0x80 };
    static const uint8_t trB[32] = {
        0x45,0x67,0x89,0xab,0xcd,0xef,0x01,0x23,0x32,0x54,0x76,0x98,0xba,0xdc,0xfe,0x10,
        0x24,0x68,0xac,0xe0,0x35,0x79,0xbd,0xf1,0x11,0x22,0x33,0x44,0x55,0x66,0x77,0x88 };
    static const uint8_t idA[] = "Alice";
    static const uint8_t idB[] = "Bob";

    sm2_affine_t G, P17, P29;
    sm2_affine_t PA_n, PB_n, PA_m, PB_m, RA_n, RA_m, RB_n, RB_m;
    uint8_t K_normal[32], K_mont[32], S1_normal[32], S2_normal[32], S1_mont[32], S2_mont[32];

    sm2_get_base_affine(&G);
    { uint8_t k17[32], k29[32];
      scalar_set_u64(k17,17); scalar_set_u64(k29,29);
      scalar_mul_normal_affine(&P17,&G,k17);
      scalar_mul_normal_affine(&P29,&G,k29); }

    printf("\n=== Comprehensive correctness suite ===\n");

    printf("=== 1) Domain conversion tests ===\n");
    if (!test_domain_conversion(&G,"G")) return 0;
    if (!test_domain_conversion(&P17,"[17]G")) return 0;
    if (!test_domain_conversion(&P29,"[29]G")) return 0;
    printf("\n");

    printf("=== 2) Curve operation diff tests ===\n");
    if (!test_curve_ops_diff(&G,&P17,"G and [17]G")) return 0;
    if (!test_curve_ops_diff(&P17,&P29,"[17]G and [29]G")) return 0;
    printf("\n");

    printf("=== 3) Scalar multiplication diff tests ===\n");
    { uint8_t k[32];
      int bits[] = {0,1,2,7,8,16,31,32,63,64,95,127,128,191,255};
      int nbits = (int)(sizeof(bits)/sizeof(bits[0]));

      scalar_set_zero(k); if (!run_scalar_diff_case("k=0, P=G",&G,k,1)) return 0;
      scalar_set_u64(k,1); if (!run_scalar_diff_case("k=1, P=G",&G,k,1)) return 0;
      scalar_set_u64(k,2); if (!run_scalar_diff_case("k=2, P=G",&G,k,1)) return 0;
      scalar_set_u64(k,3); if (!run_scalar_diff_case("k=3, P=G",&G,k,1)) return 0;
      scalar_set_u64(k,17); if (!run_scalar_diff_case("k=17, P=G",&G,k,1)) return 0;
      scalar_set_u64(k,255); if (!run_scalar_diff_case("k=255, P=G",&G,k,1)) return 0;
      for (int i=0;i<nbits;i++) {
          char name[64]; scalar_set_pow2(k,bits[i]);
          snprintf(name,sizeof(name),"k=2^%d, P=G",bits[i]);
          if (!run_scalar_diff_case(name,&G,k,1)) return 0; }
      scalar_set_u64(k,5); if (!run_scalar_diff_case("k=5, P=[17]G",&P17,k,1)) return 0;
      scalar_set_u64(k,31); if (!run_scalar_diff_case("k=31, P=[17]G",&P17,k,1)) return 0;
      for (int i=0;i<20;i++) {
          char name[64]; scalar_random_nonzero_test(k);
          snprintf(name,sizeof(name),"random scalar #%d, P=G",i+1);
          if (!run_scalar_diff_case(name,&G,k,0)) return 0; }
      for (int i=0;i<10;i++) {
          char name[64]; scalar_random_nonzero_test(k);
          snprintf(name,sizeof(name),"random scalar #%d, P=[17]G",i+1);
          if (!run_scalar_diff_case(name,&P17,k,0)) return 0; } }
    printf("\n");

    printf("=== 4) Fixed-base special tests ===\n");
    { uint8_t k[32];
      scalar_set_u64(k,1);   if (!run_fixedbase_special_case("k=1",k)) return 0;
      scalar_set_u64(k,17);  if (!run_fixedbase_special_case("k=17",k)) return 0;
      scalar_set_u64(k,255); if (!run_fixedbase_special_case("k=255",k)) return 0;
      scalar_set_pow2(k,127); if (!run_fixedbase_special_case("k=2^127",k)) return 0;
      for (int i=0;i<10;i++) {
          char name[64]; scalar_random_nonzero_test(k);
          snprintf(name,sizeof(name),"random #%d",i+1);
          if (!run_fixedbase_special_case(name,k)) return 0; } }
    printf("\n");

    printf("=== 5) KEX differential tests ===\n");
    derive_public_key_normal(&PA_n,tdA); derive_public_key_normal(&PB_n,tdB);
    derive_public_key_mont(&PA_m,tdA);   derive_public_key_mont(&PB_m,tdB);
    if (!point_equal(&PA_n,&PA_m)) { printf("FAIL: PA normal != PA mont\n"); return 0; }
    if (!point_equal(&PB_n,&PB_m)) { printf("FAIL: PB normal != PB mont\n"); return 0; }
    printf("PASS: public key derivation matches\n");

    if (!sm2_kex_initiator_gen_RA(&RA_n,trA))      { printf("FAIL: normal RA gen\n"); return 0; }
    if (!sm2_kex_initiator_gen_RA_mont(&RA_m,trA)) { printf("FAIL: mont RA gen\n");   return 0; }
    if (!point_equal(&RA_n,&RA_m)) { printf("FAIL: RA normal != RA mont\n"); return 0; }
    if (!sm2_kex_initiator_gen_RA(&RB_n,trB))      { printf("FAIL: normal RB gen\n"); return 0; }
    if (!sm2_kex_initiator_gen_RA_mont(&RB_m,trB)) { printf("FAIL: mont RB gen\n");   return 0; }
    if (!point_equal(&RB_n,&RB_m)) { printf("FAIL: RB normal != RB mont\n"); return 0; }
    printf("PASS: ephemeral public key generation matches\n");

    if (!run_one_kex_normal(K_normal,S1_normal,S2_normal,tdA,tdB,trA,trB,
                            &PA_n,&PB_n,idA,sizeof(idA)-1,idB,sizeof(idB)-1))
        { printf("FAIL: normal KEX\n"); return 0; }
    printf("PASS: normal KEX succeeded\n");

    if (!run_one_kex_mont(K_mont,S1_mont,S2_mont,tdA,tdB,trA,trB,
                          &PA_n,&PB_n,idA,sizeof(idA)-1,idB,sizeof(idB)-1))
        { printf("FAIL: Montgomery KEX\n"); return 0; }
    printf("PASS: Montgomery KEX succeeded\n");

    if (!bytes_equal(K_normal,K_mont,32))   { printf("FAIL: shared key mismatch\n"); return 0; }
    if (!bytes_equal(S1_normal,S1_mont,32)) { printf("FAIL: S1 mismatch\n"); return 0; }
    if (!bytes_equal(S2_normal,S2_mont,32)) { printf("FAIL: S2 mismatch\n"); return 0; }
    printf("PASS: normal and Montgomery outputs match\n");
    print_hex("SharedKey",K_normal,32);
    print_hex("S1",S1_normal,32);
    print_hex("S2",S2_normal,32);

    { const int rounds = 10;
      uint8_t rdA[32],rdB[32],rrA[32],rrB[32];
      sm2_affine_t RPA,RPB;
      for (int i=0;i<rounds;i++) {
          uint8_t Kn[32],Km[32],S1n[32],S2n[32],S1m[32],S2m[32];
          scalar_random_nonzero_test(rdA); scalar_random_nonzero_test(rdB);
          scalar_random_nonzero_test(rrA); scalar_random_nonzero_test(rrB);
          derive_public_key_normal(&RPA,rdA); derive_public_key_normal(&RPB,rdB);
          if (!run_one_kex_normal(Kn,S1n,S2n,rdA,rdB,rrA,rrB,&RPA,&RPB,
                                  idA,sizeof(idA)-1,idB,sizeof(idB)-1))
              { printf("FAIL: random normal KEX #%d\n",i+1); return 0; }
          if (!run_one_kex_mont(Km,S1m,S2m,rdA,rdB,rrA,rrB,&RPA,&RPB,
                                idA,sizeof(idA)-1,idB,sizeof(idB)-1))
              { printf("FAIL: random mont KEX #%d\n",i+1); return 0; }
          if (!bytes_equal(Kn,Km,32)||!bytes_equal(S1n,S1m,32)||!bytes_equal(S2n,S2m,32))
              { printf("FAIL: random KEX mismatch #%d\n",i+1); return 0; }
      }
      printf("PASS: random KEX differential tests (%d rounds)\n",rounds); }

    printf("=== Comprehensive correctness suite PASSED ===\n\n");
    return 1;
}
#endif /* USE_FP_MONT - correctness suite */

/* ============================================================
 * main
 * ============================================================ */
int main(void) {
    printf("SM2 Benchmark (strict, reproducible)\n");
    printf("=====================================\n");
    printf("iters: fp=%d, inv=%d, mul=%d, kex=%d\n\n",
           ITER_FP, ITER_INV, ITER_MUL, ITER_KEX);

    sm2_curve_init_once();

    /* 跳过正确性测试 */
    /*
    if (!run_comprehensive_correctness_suite()) {
        puts("[FAIL] comprehensive correctness suite failed");
        return 1;
    }
    */

    sm2_affine_t G, PA, PB, RA, RB;
#ifndef USE_FP_MONT
    sm2_jacobian_t J;
#endif
    uint8_t dA[32], dB[32], rA[32], rB[32];
    uint8_t KA[32], KB[32], S1[32], S2[32], SA[32], SB[32];
    const uint8_t idA[] = "Alice";
    const uint8_t idB[] = "Bob";

    sm2_get_base_affine(&G);
    fixed_scalar(dA, 0x11); fixed_scalar(dB, 0x22);
    fixed_scalar(rA, 0x33); fixed_scalar(rB, 0x44);

#ifndef USE_FP_MONT
    /* ── fp_mul ── */
    { uint64_t cyc[REPEAT_FP]; fp_t t; memset(&t,0x12,sizeof(t)); t.v[0]|=1;
      for (int r=0;r<REPEAT_FP;r++) {
          uint64_t c0=rdcycle();
          for (int i=0;i<ITER_FP;i++) sink_mem(&t,sizeof(t));
          uint64_t empty=rdcycle()-c0;
          for (int i=0;i<200;i++) fp_mul(&t,&t,&t);
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_FP;i++) { fp_mul(&t,&t,&t); if((i&127)==0) sink_mem(&t,sizeof(t)); }
          uint64_t tot=rdcycle()-s0;
          cyc[r]=(tot>empty?(tot-empty):tot)/(uint64_t)ITER_FP;
          sink_mem(&t,sizeof(t)); }
      print_stats("fp_mul (Precise)", cyc, REPEAT_FP); }

    /* ── fp_sqr (Precise) ── */
    { uint64_t cyc[REPEAT_FP]; fp_t s; memset(&s,0x56,sizeof(s)); s.v[0]|=1;
      for (int r=0;r<REPEAT_FP;r++) {
          uint64_t c0=rdcycle();
          for (int i=0;i<ITER_FP;i++) sink_mem(&s,sizeof(s));
          uint64_t empty=rdcycle()-c0;
          for (int i=0;i<200;i++) fp_sqr(&s,&s);
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_FP;i++) { fp_sqr(&s,&s); if((i&127)==0) sink_mem(&s,sizeof(s)); }
          uint64_t tot=rdcycle()-s0;
          cyc[r]=(tot>empty?(tot-empty):tot)/(uint64_t)ITER_FP;
          sink_mem(&s,sizeof(s)); }
      print_stats("fp_sqr (Precise)", cyc, REPEAT_FP); }

    /* ── fp_sqr (chained) ── */
    { uint64_t cyc[REPEAT_FP]; fp_t xs; memset(&xs,0x5a,sizeof(xs)); xs.v[0]|=1;
      for (int r=0;r<REPEAT_FP;r++) {
          for (int i=0;i<200;i++) fp_sqr(&xs,&xs);
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_FP;i++) fp_sqr(&xs,&xs);
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_FP;
          sink_mem(&xs,sizeof(xs)); }
      print_stats("fp_sqr (chained)", cyc, REPEAT_FP); }

    /* ── fp_inv ── */
    { uint64_t cyc[REPEAT_FP]; fp_t xi,ri; memset(&xi,0x9a,sizeof(xi)); xi.v[0]|=1;
      for (int r=0;r<REPEAT_FP;r++) {
          for (int i=0;i<10;i++) { fp_inv(&ri,&xi); sink_mem(&ri,sizeof(ri)); }
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_INV;i++) {
              fp_inv(&ri,&xi); fp_add(&xi,&xi,&ri); xi.v[0]|=1; }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_INV;
          sink_mem(&ri,sizeof(ri)); sink_mem(&xi,sizeof(xi)); }
      print_stats("fp_inv (avg)", cyc, REPEAT_FP); }
#endif /* !USE_FP_MONT */

#ifdef USE_FP_MONT
    /* ── fp_mont_mul ── */
    { uint64_t cyc[REPEAT_FP]; fp_t t; memset(&t,0x12,sizeof(t)); t.v[0]|=1;
      fp_to_mont(&t,&t);
      for (int r=0;r<REPEAT_FP;r++) {
          for (int i=0;i<200;i++) fp_mont_mul(&t,&t,&t);
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_FP;i++) { fp_mont_mul(&t,&t,&t); if((i&127)==0) sink_mem(&t,sizeof(t)); }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_FP;
          sink_mem(&t,sizeof(t)); }
      print_stats(FP_MONT_MUL_LABEL, cyc, REPEAT_FP); }

    /* ── fp_mont_sqr ── */
    { uint64_t cyc[REPEAT_FP]; fp_t s; memset(&s,0x56,sizeof(s)); s.v[0]|=1;
      fp_to_mont(&s,&s);
      for (int r=0;r<REPEAT_FP;r++) {
          for (int i=0;i<200;i++) fp_mont_sqr(&s,&s);
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_FP;i++) { fp_mont_sqr(&s,&s); if((i&127)==0) sink_mem(&s,sizeof(s)); }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_FP;
          sink_mem(&s,sizeof(s)); }
      print_stats(FP_MONT_SQR_LABEL, cyc, REPEAT_FP); }

    /* ── fp_mont_inv ── */
    { uint64_t cyc[REPEAT_INV]; fp_t xi,ri; memset(&xi,0x9a,sizeof(xi)); xi.v[0]|=1;
      fp_to_mont(&xi,&xi);
      for (int r=0;r<REPEAT_INV;r++) {
          for (int i=0;i<10;i++) { fp_mont_inv(&ri,&xi); sink_mem(&ri,sizeof(ri)); }
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_INV;i++) {
              fp_mont_inv(&ri,&xi); fp_mont_add(&xi,&xi,&ri); }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_INV;
          sink_mem(&ri,sizeof(ri)); sink_mem(&xi,sizeof(xi)); }
      print_stats(FP_MONT_INV_LABEL, cyc, REPEAT_INV); }
#endif /* USE_FP_MONT */

#ifdef HAVE_FN_BARRETT_CORE
#ifndef NO_TEST_FN_INTERNALS
    { uint64_t cyc[REPEAT_FN]; fn_t fa,fb; uint64_t prod[8];
      memset(&fa,0x12,sizeof(fa)); memset(&fb,0x34,sizeof(fb)); fa.v[0]|=1; fb.v[0]|=1;
      for (int r=0;r<REPEAT_FN;r++) {
          for (int i=0;i<200;i++) { fn_mul_4x4_u512_ct(prod,&fa,&fb); fa.v[0]^=prod[0]; fb.v[1]^=prod[1]; }
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_FN;i++) {
              fn_mul_4x4_u512_ct(prod,&fa,&fb);
              fa.v[0]+=prod[0]; fa.v[1]^=prod[1]; fb.v[0]+=prod[2]; fb.v[1]^=prod[3];
              if((i&255)==0) sink_mem(prod,sizeof(prod));
          }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_FN;
          sink_mem(prod,sizeof(prod)); sink_mem(&fa,sizeof(fa)); sink_mem(&fb,sizeof(fb));
      }
      print_stats("fn_mul_4x4_u512_ct (Barrett)", cyc, REPEAT_FN); }
#endif

    { uint64_t cyc[REPEAT_FN]; fn_t fm_a,fm_b,fm_r;
      memset(&fm_a,0x21,sizeof(fm_a)); memset(&fm_b,0x43,sizeof(fm_b)); fm_a.v[0]|=1; fm_b.v[0]|=1;
      for (int r=0;r<REPEAT_FN;r++) {
          for (int i=0;i<100;i++) { fn_mul(&fm_r,&fm_a,&fm_b); fm_a=fm_r; }
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_FN;i++) { fn_mul(&fm_r,&fm_a,&fm_b); fm_a=fm_r; if((i&255)==0) sink_mem(&fm_r,sizeof(fm_r)); }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_FN;
          sink_mem(&fm_r,sizeof(fm_r));
      }
      print_stats("fn_mul (Barrett API/full)", cyc, REPEAT_FN); }
#endif /* HAVE_FN_BARRETT_CORE */

#ifdef HAVE_FN_MONT_CORE
#ifndef NO_TEST_FN_INTERNALS
    { uint64_t cyc[REPEAT_FN]; fn_t fa,fb; uint64_t prod[8];
      memset(&fa,0x12,sizeof(fa)); memset(&fb,0x34,sizeof(fb)); fa.v[0]|=1; fb.v[0]|=1;
      for (int r=0;r<REPEAT_FN;r++) {
          for (int i=0;i<200;i++) { fn_mul_4x4_u512_ct(prod,&fa,&fb); fa.v[0]^=prod[0]; fb.v[1]^=prod[1]; }
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_FN;i++) { fn_mul_4x4_u512_ct(prod,&fa,&fb); fa.v[0]+=prod[0]; fa.v[1]^=prod[1]; fb.v[0]+=prod[2]; fb.v[1]^=prod[3]; }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_FN;
          sink_mem(prod,sizeof(prod)); sink_mem(&fa,sizeof(fa)); }
      print_stats("fn_mul_4x4_u512_ct", cyc, REPEAT_FN); }

    { uint64_t cyc[REPEAT_FN]; uint64_t t_red[9]; fn_t fr;
      for (int i=0;i<8;i++) t_red[i]=0x1111111111111111ULL*(uint64_t)(i+1); t_red[8]=0;
      for (int r=0;r<REPEAT_FN;r++) {
          for (int i=0;i<200;i++) { fn_mont_reduce_ct(&fr,t_red); t_red[0]^=fr.v[0]; t_red[1]+=fr.v[1]; t_red[8]=0; }
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_FN;i++) { fn_mont_reduce_ct(&fr,t_red); t_red[0]+=fr.v[0]; t_red[1]^=fr.v[1]; t_red[2]+=fr.v[2]; t_red[3]^=fr.v[3]; t_red[8]=0; }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_FN;
          sink_mem(&fr,sizeof(fr)); }
      print_stats("fn_mont_reduce_ct", cyc, REPEAT_FN); }
#endif

    { uint64_t cyc[REPEAT_FN]; fn_t fm_a,fm_b,fm_r;
      memset(&fm_a,0x21,sizeof(fm_a)); memset(&fm_b,0x43,sizeof(fm_b)); fm_a.v[0]|=1; fm_b.v[0]|=1;
      for (int r=0;r<REPEAT_FN;r++) {
          for (int i=0;i<100;i++) { fn_mont_mul_ct(&fm_r,&fm_a,&fm_b); fm_a=fm_r; }
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_FN;i++) { fn_mont_mul_ct(&fm_r,&fm_a,&fm_b); fm_a=fm_r; if((i&255)==0) sink_mem(&fm_r,sizeof(fm_r)); }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_FN;
          sink_mem(&fm_r,sizeof(fm_r)); }
      print_stats("fn_modmul_raw (Mont)", cyc, REPEAT_FN); }

    { uint64_t cyc[REPEAT_FN]; fn_t fm_a,fm_b,fm_r;
      memset(&fm_a,0x21,sizeof(fm_a)); memset(&fm_b,0x43,sizeof(fm_b)); fm_a.v[0]|=1; fm_b.v[0]|=1;
      for (int r=0;r<REPEAT_FN;r++) {
          for (int i=0;i<100;i++) { fn_mul(&fm_r,&fm_a,&fm_b); fm_a=fm_r; }
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_FN;i++) { fn_mul(&fm_r,&fm_a,&fm_b); fm_a=fm_r; if((i&255)==0) sink_mem(&fm_r,sizeof(fm_r)); }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_FN;
          sink_mem(&fm_r,sizeof(fm_r)); }
      print_stats("fn_mul (Mont API/full)", cyc, REPEAT_FN); }
#endif /* HAVE_FN_MONT_CORE */

#ifndef USE_FP_MONT
    /* ── ScalarMul SchemeB (normal) ── */
    { uint64_t cyc[REPEAT_MUL]; uint8_t d[32];
      for (int r=0;r<REPEAT_MUL;r++) {
          rand_scalar_mod_n(d);
          for (int i=0;i<10;i++) sm2_scalar_mul_window_ct_schemeB(&J,&G,d);
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_MUL;i++) { sm2_scalar_mul_window_ct_schemeB(&J,&G,d); sink_mem(&J,sizeof(J)); }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_MUL; }
      print_stats("ScalarMul SchemeB normal ([d]G)", cyc, REPEAT_MUL); }
#endif /* !USE_FP_MONT */

#ifdef USE_FP_MONT
    /* ── ScalarMul SchemeB (Montgomery) ── */
    { uint64_t cyc[REPEAT_MUL]; uint8_t d[32];
      sm2_affine_t Gm;
      sm2_jacobian_t RmJ;
      sm2_affine_to_mont(&Gm, &G);
      for (int r=0;r<REPEAT_MUL;r++) {
          rand_scalar_mod_n(d);
          for (int i=0;i<10;i++) sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&RmJ,&Gm,d);
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_MUL;i++) { sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&RmJ,&Gm,d); sink_mem(&RmJ,sizeof(RmJ)); }
          cyc[r]=(rdcycle()-s0)/(uint64_t)ITER_MUL;
      }
      print_stats("ScalarMul SchemeB Mont-core ([d]G)", cyc, REPEAT_MUL); }
#endif /* USE_FP_MONT */

#ifndef USE_FP_MONT
    /* ── KEX normal ── */
    { uint64_t kex_gen[REPEAT_KEX], kex_resp[REPEAT_KEX], kex_init[REPEAT_KEX];
      for (int r=0;r<REPEAT_KEX;r++) {
          rand_scalar_mod_n(dA); rand_scalar_mod_n(dB);
          rand_scalar_mod_n(rA); rand_scalar_mod_n(rB);
          sm2_scalar_mul_window_ct_schemeB(&J,&G,dA); sm2_jacobian_to_affine(&PA,&J);
          sm2_scalar_mul_window_ct_schemeB(&J,&G,dB); sm2_jacobian_to_affine(&PB,&J);

          for (int i=0;i<5;i++) sm2_kex_initiator_gen_RA(&RA,rA);
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_KEX;i++) { sm2_kex_initiator_gen_RA(&RA,rA); sink_mem(&RA,sizeof(RA)); }
          kex_gen[r]=(rdcycle()-s0)/(uint64_t)ITER_KEX;

          if (!sm2_kex_initiator_gen_RA(&RA,rA)||!sm2_kex_initiator_gen_RA(&RB,rB)) return 1;
          if (!sm2_kex_initiator_compute_key(KA,32,S1,SA,NULL,dA,&PA,idA,sizeof(idA)-1,rA,&RA,&PB,idB,sizeof(idB)-1,&RB)) return 1;
          for (int i=0;i<3;i++) sm2_kex_responder_compute_key(&RB,KB,32,S2,SB,S1,dB,&PB,idB,sizeof(idB)-1,rB,&RA,&PA,idA,sizeof(idA)-1);
          s0=rdcycle();
          for (int i=0;i<ITER_KEX;i++) { sm2_kex_responder_compute_key(&RB,KB,32,S2,SB,S1,dB,&PB,idB,sizeof(idB)-1,rB,&RA,&PA,idA,sizeof(idA)-1); sink_mem(KB,32); }
          kex_resp[r]=(rdcycle()-s0)/(uint64_t)ITER_KEX;

          for (int i=0;i<3;i++) sm2_kex_initiator_compute_key(KA,32,S1,SA,S2,dA,&PA,idA,sizeof(idA)-1,rA,&RA,&PB,idB,sizeof(idB)-1,&RB);
          s0=rdcycle();
          for (int i=0;i<ITER_KEX;i++) { sm2_kex_initiator_compute_key(KA,32,S1,SA,S2,dA,&PA,idA,sizeof(idA)-1,rA,&RA,&PB,idB,sizeof(idB)-1,&RB); sink_mem(KA,32); }
          kex_init[r]=(rdcycle()-s0)/(uint64_t)ITER_KEX;
          sink_mem(S1,32); sink_mem(S2,32); }
      print_stats("KEX Initiator Gen RA (normal)",   kex_gen,  REPEAT_KEX);
      print_stats("KEX Responder normal (full)",      kex_resp, REPEAT_KEX);
      print_stats("KEX Initiator normal (verify)",    kex_init, REPEAT_KEX); }
#endif /* !USE_FP_MONT */

#ifdef USE_FP_MONT
    /* ── KEX Montgomery ── */
    { uint64_t kex_gen[REPEAT_KEX], kex_resp[REPEAT_KEX], kex_init[REPEAT_KEX];
      for (int r=0;r<REPEAT_KEX;r++) {
          rand_scalar_mod_n(dA); rand_scalar_mod_n(dB);
          rand_scalar_mod_n(rA); rand_scalar_mod_n(rB);
          sm2_scalar_mul_window_ct_schemeB_mont_to_affine(&PA,&G,dA);
          sm2_scalar_mul_window_ct_schemeB_mont_to_affine(&PB,&G,dB);

          for (int i=0;i<5;i++) sm2_kex_initiator_gen_RA_mont(&RA,rA);
          uint64_t s0=rdcycle();
          for (int i=0;i<ITER_KEX;i++) { sm2_kex_initiator_gen_RA_mont(&RA,rA); sink_mem(&RA,sizeof(RA)); }
          kex_gen[r]=(rdcycle()-s0)/(uint64_t)ITER_KEX;

          if (!sm2_kex_initiator_gen_RA_mont(&RA,rA)||!sm2_kex_initiator_gen_RA_mont(&RB,rB)) return 1;
          if (!sm2_kex_initiator_compute_key_mont(KA,32,S1,SA,NULL,dA,&PA,idA,sizeof(idA)-1,rA,&RA,&PB,idB,sizeof(idB)-1,&RB)) return 1;
          for (int i=0;i<3;i++) sm2_kex_responder_compute_key_mont(&RB,KB,32,S2,SB,S1,dB,&PB,idB,sizeof(idB)-1,rB,&RA,&PA,idA,sizeof(idA)-1);
          s0=rdcycle();
          for (int i=0;i<ITER_KEX;i++) { sm2_kex_responder_compute_key_mont(&RB,KB,32,S2,SB,S1,dB,&PB,idB,sizeof(idB)-1,rB,&RA,&PA,idA,sizeof(idA)-1); sink_mem(KB,32); }
          kex_resp[r]=(rdcycle()-s0)/(uint64_t)ITER_KEX;

          for (int i=0;i<3;i++) sm2_kex_initiator_compute_key_mont(KA,32,S1,SA,S2,dA,&PA,idA,sizeof(idA)-1,rA,&RA,&PB,idB,sizeof(idB)-1,&RB);
          s0=rdcycle();
          for (int i=0;i<ITER_KEX;i++) { sm2_kex_initiator_compute_key_mont(KA,32,S1,SA,S2,dA,&PA,idA,sizeof(idA)-1,rA,&RA,&PB,idB,sizeof(idB)-1,&RB); sink_mem(KA,32); }
          kex_init[r]=(rdcycle()-s0)/(uint64_t)ITER_KEX;
          sink_mem(S1,32); sink_mem(S2,32); }
      print_stats("KEX Initiator Gen RA (Mont-C)",   kex_gen,  REPEAT_KEX);
      print_stats("KEX Responder Mont-C (full)",      kex_resp, REPEAT_KEX);
      print_stats("KEX Initiator Mont-C (verify)",    kex_init, REPEAT_KEX); }
#endif /* USE_FP_MONT */

    printf("\n[sink] %u\n", bench_sink);
    return 0;
}
