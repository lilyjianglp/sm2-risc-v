#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

/* self */
#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/fp.h"
#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/fn.h"

#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/sm2_curve.h"
#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/sm2_scalar.h"

/* gmssl */
#include "/home/sm2/GmSSL/include/gmssl/sm2.h"
#include "/home/sm2/GmSSL/include/gmssl/sm2_z256.h"

#define SINK_STRIDE 128

#define REPEAT_FP   5
#define REPEAT_MUL  5
#define REPEAT_KEX  5
#define REPEAT_FN   5

#define ITER_FP     100000
#define ITER_MUL    1000
#define ITER_KEX    200
#define ITER_INV    200
#define ITER_FN     20000

typedef struct {
    uint64_t cyc_mean, cyc_min, cyc_max;
    uint64_t ns_mean, ns_min, ns_max;
} bench_result_t;

typedef struct {
    double cyc_mean, cyc_min, cyc_max;
    double ns_mean, ns_min, ns_max;
} bench_precise_result_t;

static bench_result_t zero = {0};
volatile uint8_t bench_sink = 0;

/* ========= prefixed self symbols =========
 * These names come from objcopy --prefix-symbols=c_ / asm_
 * Only declare the minimal set we use here.
 */
void c_fp_to_mont(fp_t *r, const fp_t *a);
void c_fp_from_mont(fp_t *r, const fp_t *a);
void c_fp_mont_mul(fp_t *r, const fp_t *a, const fp_t *b);
void c_fp_mont_sqr(fp_t *r, const fp_t *a);
void c_fp_mont_add(fp_t *r, const fp_t *a, const fp_t *b);
void c_fp_mont_sub(fp_t *r, const fp_t *a, const fp_t *b);
void c_fp_mont_inv(fp_t *r, const fp_t *a);

void asm_fp_to_mont(fp_t *r, const fp_t *a);
void asm_fp_from_mont(fp_t *r, const fp_t *a);
void asm_fp_mont_mul(fp_t *r, const fp_t *a, const fp_t *b);
void asm_fp_mont_sqr(fp_t *r, const fp_t *a);
void asm_fp_mont_add(fp_t *r, const fp_t *a, const fp_t *b);
void asm_fp_mont_sub(fp_t *r, const fp_t *a, const fp_t *b);
void asm_fp_mont_inv(fp_t *r, const fp_t *a);


void c_fn_to_mont(fn_t *r, const fn_t *a);
void c_fn_mont_mul_ct(fn_t *r, const fn_t *a, const fn_t *b);
void c_fn_mul(fn_t *r, const fn_t *a, const fn_t *b);

void asm_fn_to_mont(fn_t *r, const fn_t *a);
void asm_fn_mont_mul_ct(fn_t *r, const fn_t *a, const fn_t *b);
void asm_fn_mul(fn_t *r, const fn_t *a, const fn_t *b);

void c_sm2_get_base_affine_mont(sm2_affine_t *G);
void asm_sm2_get_base_affine_mont(sm2_affine_t *G);

void c_sm2_scalar_mul_window_ct_schemeB_mont_to_affine(
    sm2_affine_t *R, const sm2_affine_t *P, const uint8_t k[32]);

void asm_sm2_scalar_mul_window_ct_schemeB_mont_to_affine(
    sm2_affine_t *R, const sm2_affine_t *P, const uint8_t k[32]);

void c_sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(
    sm2_jacobian_t *R, const sm2_affine_t *P, const uint8_t k[32]);

void asm_sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(
    sm2_jacobian_t *R, const sm2_affine_t *P, const uint8_t k[32]);

/* ========= utils ========= */

static inline void sink_mem(const void *p, size_t n) {
    const volatile uint8_t *b = (const volatile uint8_t *)p;
    for (size_t i = 0; i < n; i++) {
        bench_sink ^= b[i];
    }
}

static inline void sink_sparse(const void *p, size_t n, int i) {
    if ((i & (SINK_STRIDE - 1)) == 0) {
        sink_mem(p, n);
    }
}

static inline uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

static inline uint64_t rdcycle(void) {
    uint64_t c;
    asm volatile ("rdcycle %0" : "=r"(c));
    return c;
}

static inline uint64_t div_round_down(uint64_t x, uint64_t y) {
    return y ? (x / y) : 0;
}

static bench_result_t make_result(const uint64_t *cyc, const uint64_t *ns, int n) {
    bench_result_t r;
    uint64_t cyc_sum = cyc[0], cyc_min = cyc[0], cyc_max = cyc[0];
    uint64_t ns_sum  = ns[0],  ns_min  = ns[0],  ns_max  = ns[0];

    for (int i = 1; i < n; i++) {
        cyc_sum += cyc[i];
        ns_sum  += ns[i];

        if (cyc[i] < cyc_min) cyc_min = cyc[i];
        if (cyc[i] > cyc_max) cyc_max = cyc[i];
        if (ns[i]  < ns_min)  ns_min  = ns[i];
        if (ns[i]  > ns_max)  ns_max  = ns[i];
    }

    r.cyc_mean = cyc_sum / (uint64_t)n;
    r.cyc_min  = cyc_min;
    r.cyc_max  = cyc_max;
    r.ns_mean  = ns_sum / (uint64_t)n;
    r.ns_min   = ns_min;
    r.ns_max   = ns_max;
    return r;
}

static bench_precise_result_t make_precise_result(const double *cyc, const double *ns, int n) {
    bench_precise_result_t r;
    double cyc_sum = cyc[0], cyc_min = cyc[0], cyc_max = cyc[0];
    double ns_sum  = ns[0],  ns_min  = ns[0],  ns_max  = ns[0];

    for (int i = 1; i < n; i++) {
        cyc_sum += cyc[i];
        ns_sum  += ns[i];

        if (cyc[i] < cyc_min) cyc_min = cyc[i];
        if (cyc[i] > cyc_max) cyc_max = cyc[i];
        if (ns[i]  < ns_min)  ns_min  = ns[i];
        if (ns[i]  > ns_max)  ns_max  = ns[i];
    }

    r.cyc_mean = cyc_sum / (double)n;
    r.cyc_min  = cyc_min;
    r.cyc_max  = cyc_max;
    r.ns_mean  = ns_sum / (double)n;
    r.ns_min   = ns_min;
    r.ns_max   = ns_max;
    return r;
}

static void print_row_3way(const char *name,
                           bench_result_t c_r,
                           bench_result_t asm_r,
                           bench_result_t gmssl_r) {
    double asm_vs_c     = (double)c_r.cyc_mean / (double)asm_r.cyc_mean;
    double c_vs_gmssl   = (double)gmssl_r.cyc_mean > 0 ? (double)gmssl_r.cyc_mean / (double)c_r.cyc_mean : 0.0;
    double asm_vs_gmssl = (double)gmssl_r.cyc_mean > 0 ? (double)gmssl_r.cyc_mean / (double)asm_r.cyc_mean : 0.0;

    printf("%-16s  c=%-8llu  asm=%-8llu  gmssl=%-8llu  "
           "[speedup] asm>c %.2fx  c>gmssl %.2fx  asm>gmssl %.2fx\n",
           name,
           (unsigned long long)c_r.cyc_mean,
           (unsigned long long)asm_r.cyc_mean,
           (unsigned long long)gmssl_r.cyc_mean,
           asm_vs_c,
           c_vs_gmssl,
           asm_vs_gmssl);
}

static void print_row_3way_precise(const char *name,
                                   bench_precise_result_t c_r,
                                   bench_precise_result_t asm_r,
                                   bench_precise_result_t gmssl_r) {
    double asm_vs_c     = c_r.cyc_mean / asm_r.cyc_mean;
    double c_vs_gmssl   = gmssl_r.cyc_mean > 0.0 ? gmssl_r.cyc_mean / c_r.cyc_mean : 0.0;
    double asm_vs_gmssl = gmssl_r.cyc_mean > 0.0 ? gmssl_r.cyc_mean / asm_r.cyc_mean : 0.0;

    printf("%-16s  c=%-8.2f  asm=%-8.2f  gmssl=%-8.2f  "
           "[speedup] asm>c %.3fx  c>gmssl %.3fx  asm>gmssl %.3fx\n",
           name,
           c_r.cyc_mean,
           asm_r.cyc_mean,
           gmssl_r.cyc_mean,
           asm_vs_c,
           c_vs_gmssl,
           asm_vs_gmssl);
    printf("%-16s  c[min,max]=[%.2f, %.2f]  asm[min,max]=[%.2f, %.2f]  gmssl[min,max]=[%.2f, %.2f]\n",
           "",
           c_r.cyc_min, c_r.cyc_max,
           asm_r.cyc_min, asm_r.cyc_max,
           gmssl_r.cyc_min, gmssl_r.cyc_max);
}
static void fixed_scalar(uint8_t k[32], uint8_t v) {
    memset(k, v, 32);
    k[31] |= 1;
}

static void z256_from_be32(sm2_z256_t out, const uint8_t in[32]) {
    for (int i = 0; i < 4; i++) {
        uint64_t w = 0;
        for (int j = 0; j < 8; j++) {
            w = (w << 8) | in[i * 8 + j];
        }
        out[3 - i] = w;
    }
}

/* ========= self-c fp_mont_mul ========= */

static bench_result_t bench_self_c_fp_mont_mul(void) {
    fp_t a, b, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    c_fp_to_mont(&a, &a);
    c_fp_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            c_fp_mont_mul(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            c_fp_mont_mul(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= self-asm fp_mont_mul ========= */

static bench_result_t bench_self_asm_fp_mont_mul(void) {
    fp_t a, b, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    asm_fp_to_mont(&a, &a);
    asm_fp_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            asm_fp_mont_mul(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            asm_fp_mont_mul(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= gmssl fp_mont_mul ========= */

static bench_result_t bench_gmssl_modp_mont_mul(void) {
    sm2_z256_t a, b, r;
    uint8_t a_bytes[32], b_bytes[32];
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    fixed_scalar(a_bytes, 0x12);
    fixed_scalar(b_bytes, 0x34);
    z256_from_be32(a, a_bytes);
    z256_from_be32(b, b_bytes);

    sm2_z256_modp_to_mont(a, a);
    sm2_z256_modp_to_mont(b, b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        sm2_z256_t aa, bb, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        memcpy(aa, a, sizeof(aa));
        memcpy(bb, b, sizeof(bb));

        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_mont_mul(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            memcpy(tmp, aa, sizeof(tmp));
            sink_sparse(tmp, sizeof(tmp), i);
            memcpy(aa, tmp, sizeof(tmp));
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            sm2_z256_modp_mont_mul(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= self-c fp_mont_sqr ========= */

static bench_result_t bench_self_c_fp_mont_sqr(void) {
    fp_t a, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;

    c_fp_to_mont(&a, &a);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            c_fp_mont_sqr(&r, &aa);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            c_fp_mont_sqr(&r, &aa);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= self-asm fp_mont_sqr ========= */

static bench_result_t bench_self_asm_fp_mont_sqr(void) {
    fp_t a, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;

    asm_fp_to_mont(&a, &a);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            asm_fp_mont_sqr(&r, &aa);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            asm_fp_mont_sqr(&r, &aa);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= gmssl fp_mont_sqr ========= */

static bench_result_t bench_gmssl_modp_mont_sqr(void) {
    sm2_z256_t a, r;
    uint8_t a_bytes[32];
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    fixed_scalar(a_bytes, 0x12);
    z256_from_be32(a, a_bytes);

    sm2_z256_modp_to_mont(a, a);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        sm2_z256_t aa, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        memcpy(aa, a, sizeof(aa));

        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_mont_sqr(r, aa);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            memcpy(tmp, aa, sizeof(tmp));
            sink_sparse(tmp, sizeof(tmp), i);
            memcpy(aa, tmp, sizeof(tmp));
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            sm2_z256_modp_mont_sqr(r, aa);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}
/* ========= self-c fp_mont_add ========= */

static bench_result_t bench_self_c_fp_mont_add(void) {
    fp_t a, b, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    c_fp_to_mont(&a, &a);
    c_fp_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            c_fp_mont_add(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            c_fp_mont_add(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= self-asm fp_mont_add ========= */

static bench_result_t bench_self_asm_fp_mont_add(void) {
    fp_t a, b, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    asm_fp_to_mont(&a, &a);
    asm_fp_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            asm_fp_mont_add(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            asm_fp_mont_add(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

static bench_precise_result_t bench_self_c_fp_mont_add_precise(void) {
    fp_t a, b, r;
    double cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    c_fp_to_mont(&a, &a);
    c_fp_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;
        uint64_t net_cyc, net_ns;

        for (int i = 0; i < 200; i++) {
            c_fp_mont_add(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            c_fp_mont_add(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        net_cyc = full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc;
        net_ns  = full_ns > empty_ns ? (full_ns - empty_ns) : full_ns;

        cyc[rr] = (double)net_cyc / (double)ITER_FP;
        ns[rr]  = (double)net_ns  / (double)ITER_FP;

        sink_mem(&aa, sizeof(aa));
    }

    return make_precise_result(cyc, ns, REPEAT_FP);
}

static bench_precise_result_t bench_self_asm_fp_mont_add_precise(void) {
    fp_t a, b, r;
    double cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    asm_fp_to_mont(&a, &a);
    asm_fp_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;
        uint64_t net_cyc, net_ns;

        for (int i = 0; i < 200; i++) {
            asm_fp_mont_add(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            asm_fp_mont_add(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        net_cyc = full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc;
        net_ns  = full_ns > empty_ns ? (full_ns - empty_ns) : full_ns;

        cyc[rr] = (double)net_cyc / (double)ITER_FP;
        ns[rr]  = (double)net_ns  / (double)ITER_FP;

        sink_mem(&aa, sizeof(aa));
    }

    return make_precise_result(cyc, ns, REPEAT_FP);
}

/* ========= gmssl fp_mont_add ========= */

static bench_result_t bench_gmssl_modp_add(void) {
    sm2_z256_t a, b, r;
    uint8_t a_bytes[32], b_bytes[32];
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    fixed_scalar(a_bytes, 0x12);
    fixed_scalar(b_bytes, 0x34);
    z256_from_be32(a, a_bytes);
    z256_from_be32(b, b_bytes);

    sm2_z256_modp_to_mont(a, a);
    sm2_z256_modp_to_mont(b, b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        sm2_z256_t aa, bb, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        memcpy(aa, a, sizeof(aa));
        memcpy(bb, b, sizeof(bb));

        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_add(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            memcpy(tmp, aa, sizeof(tmp));
            sink_sparse(tmp, sizeof(tmp), i);
            memcpy(aa, tmp, sizeof(tmp));
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            sm2_z256_modp_add(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

static bench_precise_result_t bench_gmssl_modp_add_precise(void) {
    sm2_z256_t a, b, r;
    uint8_t a_bytes[32], b_bytes[32];
    double cyc[REPEAT_FP], ns[REPEAT_FP];

    fixed_scalar(a_bytes, 0x12);
    fixed_scalar(b_bytes, 0x34);
    z256_from_be32(a, a_bytes);
    z256_from_be32(b, b_bytes);

    sm2_z256_modp_to_mont(a, a);
    sm2_z256_modp_to_mont(b, b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        sm2_z256_t aa, bb, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;
        uint64_t net_cyc, net_ns;

        memcpy(aa, a, sizeof(aa));
        memcpy(bb, b, sizeof(bb));

        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_add(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            memcpy(tmp, aa, sizeof(tmp));
            sink_sparse(tmp, sizeof(tmp), i);
            memcpy(aa, tmp, sizeof(tmp));
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            sm2_z256_modp_add(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        net_cyc = full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc;
        net_ns  = full_ns > empty_ns ? (full_ns - empty_ns) : full_ns;

        cyc[rr] = (double)net_cyc / (double)ITER_FP;
        ns[rr]  = (double)net_ns  / (double)ITER_FP;

        sink_mem(aa, sizeof(aa));
    }

    return make_precise_result(cyc, ns, REPEAT_FP);
}

/* ========= self-c fp_mont_sub ========= */

static bench_result_t bench_self_c_fp_mont_sub(void) {
    fp_t a, b, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x56, sizeof(a));
    memset(&b, 0x12, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    c_fp_to_mont(&a, &a);
    c_fp_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            c_fp_mont_sub(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            c_fp_mont_sub(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= self-asm fp_mont_sub ========= */

static bench_result_t bench_self_asm_fp_mont_sub(void) {
    fp_t a, b, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x56, sizeof(a));
    memset(&b, 0x12, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    asm_fp_to_mont(&a, &a);
    asm_fp_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            asm_fp_mont_sub(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            asm_fp_mont_sub(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= gmssl fp_mont_sub ========= */

static bench_result_t bench_gmssl_modp_sub(void) {
    sm2_z256_t a, b, r;
    uint8_t a_bytes[32], b_bytes[32];
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    fixed_scalar(a_bytes, 0x12);
    fixed_scalar(b_bytes, 0x34);
    z256_from_be32(a, a_bytes);
    z256_from_be32(b, b_bytes);

    sm2_z256_modp_to_mont(a, a);
    sm2_z256_modp_to_mont(b, b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        sm2_z256_t aa, bb, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        memcpy(aa, a, sizeof(aa));
        memcpy(bb, b, sizeof(bb));

        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_sub(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            memcpy(tmp, aa, sizeof(tmp));
            sink_sparse(tmp, sizeof(tmp), i);
            memcpy(aa, tmp, sizeof(tmp));
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            sm2_z256_modp_sub(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= self-c fp_to_mont ========= */

static bench_result_t bench_self_c_fp_to_mont(void) {
    fp_t a, r, tmp;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            c_fp_to_mont(&r, &a);
            sink_sparse(&r, sizeof(r), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = a;
            sink_sparse(&tmp, sizeof(tmp), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            c_fp_to_mont(&r, &a);
            sink_sparse(&r, sizeof(r), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&r, sizeof(r));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= self-asm fp_to_mont ========= */

static bench_result_t bench_self_asm_fp_to_mont(void) {
    fp_t a, r, tmp;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            asm_fp_to_mont(&r, &a);
            sink_sparse(&r, sizeof(r), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = a;
            sink_sparse(&tmp, sizeof(tmp), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            asm_fp_to_mont(&r, &a);
            sink_sparse(&r, sizeof(r), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&r, sizeof(r));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= gmssl fp_to_mont ========= */

static bench_result_t bench_gmssl_modp_to_mont(void) {
    sm2_z256_t a, r, tmp;
    uint8_t a_bytes[32];
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    fixed_scalar(a_bytes, 0x12);
    z256_from_be32(a, a_bytes);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_to_mont(a, r);
            sink_sparse(r, sizeof(r), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            memcpy(tmp, a, sizeof(tmp));
            sink_sparse(tmp, sizeof(tmp), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            sm2_z256_modp_to_mont(a, r);
            sink_sparse(r, sizeof(r), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(r, sizeof(r));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= self-c fp_from_mont ========= */

static bench_result_t bench_self_c_fp_from_mont(void) {
    fp_t a, r, tmp;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;
    c_fp_to_mont(&a, &a);   /* 固定 Montgomery 域输入 */

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            c_fp_from_mont(&r, &a);
            sink_sparse(&r, sizeof(r), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = a;
            sink_sparse(&tmp, sizeof(tmp), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            c_fp_from_mont(&r, &a);
            sink_sparse(&r, sizeof(r), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&r, sizeof(r));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= self-asm fp_from_mont ========= */

static bench_result_t bench_self_asm_fp_from_mont(void) {
    fp_t a, r, tmp;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;
    asm_fp_to_mont(&a, &a);   /* 固定 Montgomery 域输入 */

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            asm_fp_from_mont(&r, &a);
            sink_sparse(&r, sizeof(r), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            tmp = a;
            sink_sparse(&tmp, sizeof(tmp), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            asm_fp_from_mont(&r, &a);
            sink_sparse(&r, sizeof(r), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(&r, sizeof(r));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= gmssl fp_from_mont ========= */

static bench_result_t bench_gmssl_modp_from_mont(void) {
    sm2_z256_t a, r, tmp;
    uint8_t a_bytes[32];
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    fixed_scalar(a_bytes, 0x12);
    z256_from_be32(a, a_bytes);
    sm2_z256_modp_to_mont(a, a);   /* 固定 Montgomery 域输入 */

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_from_mont(r, a);
            sink_sparse(r, sizeof(r), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            memcpy(tmp, a, sizeof(tmp));
            sink_sparse(tmp, sizeof(tmp), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            sm2_z256_modp_from_mont(r, a);
            sink_sparse(r, sizeof(r), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FP
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FP
        );

        sink_mem(r, sizeof(r));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= self-c fp_mont_inv ========= */

static bench_result_t bench_self_c_fp_mont_inv(void) {
    fp_t a, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;

    c_fp_to_mont(&a, &a);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 20; i++) {
            c_fp_mont_inv(&r, &aa);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_INV; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_INV; i++) {
            c_fp_mont_inv(&r, &aa);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_INV
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_INV
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}
/* ========= self-asm fp_mont_inv ========= */

static bench_result_t bench_self_asm_fp_mont_inv(void) {
    fp_t a, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;

    asm_fp_to_mont(&a, &a);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 20; i++) {
            asm_fp_mont_inv(&r, &aa);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_INV; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_INV; i++) {
            asm_fp_mont_inv(&r, &aa);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_INV
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_INV
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}



/* ========= gmssl fp_mont_inv ========= */

static bench_result_t bench_gmssl_modp_mont_inv(void) {
    sm2_z256_t a, r;
    uint8_t a_bytes[32];
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    fixed_scalar(a_bytes, 0x12);
    z256_from_be32(a, a_bytes);

    sm2_z256_modp_to_mont(a, a);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        sm2_z256_t aa, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        memcpy(aa, a, sizeof(aa));

        for (int i = 0; i < 20; i++) {
            sm2_z256_modp_mont_inv(r, aa);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_INV; i++) {
            memcpy(tmp, aa, sizeof(tmp));
            sink_sparse(tmp, sizeof(tmp), i);
            memcpy(aa, tmp, sizeof(tmp));
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_INV; i++) {
            sm2_z256_modp_mont_inv(r, aa);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_INV
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_INV
        );

        sink_mem(aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FP);
}

/* ========= self-c fn_mont_mul ========= */

static bench_result_t bench_self_c_fn_mont_mul(void) {
    fn_t a, b, r;
    uint64_t cyc[REPEAT_FN], ns[REPEAT_FN];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    c_fn_to_mont(&a, &a);
    c_fn_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FN; rr++) {
        fn_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            c_fn_mont_mul_ct(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            c_fn_mont_mul_ct(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FN
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FN
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FN);
}

/* ========= self-asm fn_mont_mul ========= */

static bench_result_t bench_self_asm_fn_mont_mul(void) {
    fn_t a, b, r;
    uint64_t cyc[REPEAT_FN], ns[REPEAT_FN];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    asm_fn_to_mont(&a, &a);
    asm_fn_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FN; rr++) {
        fn_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            asm_fn_mont_mul_ct(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            asm_fn_mont_mul_ct(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FN
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FN
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FN);
}

/* ========= gmssl fn_mont_mul ========= */

static bench_result_t bench_gmssl_modn_mont_mul(void) {
    sm2_z256_t a, b, r;
    uint8_t a_bytes[32], b_bytes[32];
    uint64_t cyc[REPEAT_FN], ns[REPEAT_FN];

    fixed_scalar(a_bytes, 0x12);
    fixed_scalar(b_bytes, 0x34);
    z256_from_be32(a, a_bytes);
    z256_from_be32(b, b_bytes);

    sm2_z256_modn_to_mont(a, a);
    sm2_z256_modn_to_mont(b, b);

    for (int rr = 0; rr < REPEAT_FN; rr++) {
        sm2_z256_t aa, bb, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        memcpy(aa, a, sizeof(aa));
        memcpy(bb, b, sizeof(bb));

        for (int i = 0; i < 200; i++) {
            sm2_z256_modn_mont_mul(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            memcpy(tmp, aa, sizeof(tmp));
            sink_sparse(tmp, sizeof(tmp), i);
            memcpy(aa, tmp, sizeof(tmp));
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            sm2_z256_modn_mont_mul(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FN
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FN
        );

        sink_mem(aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FN);
}

/* ========= self-c fn_mul ========= */

static bench_result_t bench_self_c_fn_mul(void) {
    fn_t a, b, r;
    uint64_t cyc[REPEAT_FN], ns[REPEAT_FN];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    for (int rr = 0; rr < REPEAT_FN; rr++) {
        fn_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            c_fn_mul(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            c_fn_mul(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FN
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FN
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FN);
}
/* ========= self-asm fn_mul ========= */

static bench_result_t bench_self_asm_fn_mul(void) {
    fn_t a, b, r;
    uint64_t cyc[REPEAT_FN], ns[REPEAT_FN];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    for (int rr = 0; rr < REPEAT_FN; rr++) {
        fn_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 200; i++) {
            asm_fn_mul(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            tmp = aa;
            sink_sparse(&tmp, sizeof(tmp), i);
            aa = tmp;
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            asm_fn_mul(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FN
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FN
        );

        sink_mem(&aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FN);
}

/* ========= gmssl fn_mul ========= */

static bench_result_t bench_gmssl_modn_mul(void) {
    sm2_z256_t a, b, r;
    uint8_t a_bytes[32], b_bytes[32];
    uint64_t cyc[REPEAT_FN], ns[REPEAT_FN];

    fixed_scalar(a_bytes, 0x12);
    fixed_scalar(b_bytes, 0x34);
    z256_from_be32(a, a_bytes);
    z256_from_be32(b, b_bytes);

    sm2_z256_modn_to_mont(a, a);
    sm2_z256_modn_to_mont(b, b);

    for (int rr = 0; rr < REPEAT_FN; rr++) {
        sm2_z256_t aa, bb, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        memcpy(aa, a, sizeof(aa));
        memcpy(bb, b, sizeof(bb));

        for (int i = 0; i < 200; i++) {
            sm2_z256_modn_mul(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            memcpy(tmp, aa, sizeof(tmp));
            sink_sparse(tmp, sizeof(tmp), i);
            memcpy(aa, tmp, sizeof(tmp));
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            sm2_z256_modn_mul(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_FN
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_FN
        );

        sink_mem(aa, sizeof(aa));
    }

    return make_result(cyc, ns, REPEAT_FN);
}
/* ========= self-c point_mul_G ========= */

static bench_result_t bench_self_c_point_mul_generator(void) {
    sm2_affine_t Gm;
    sm2_jacobian_t R, T;
    uint8_t d[32];
    uint64_t cyc[REPEAT_MUL], ns[REPEAT_MUL];

    c_sm2_get_base_affine_mont(&Gm);
    fixed_scalar(d, 0x12);

    for (int rr = 0; rr < REPEAT_MUL; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 20; i++) {
            c_sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Gm, d);
            sink_sparse(&R, sizeof(R), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            T = R;
            sink_sparse(&T, sizeof(T), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            c_sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Gm, d);
            sink_sparse(&R, sizeof(R), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_MUL
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_MUL
        );

        sink_mem(&R, sizeof(R));
    }

    return make_result(cyc, ns, REPEAT_MUL);
}

/* ========= self-asm point_mul_G ========= */

static bench_result_t bench_self_asm_point_mul_generator(void) {
    sm2_affine_t Gm;
    sm2_jacobian_t R, T;
    uint8_t d[32];
    uint64_t cyc[REPEAT_MUL], ns[REPEAT_MUL];

    asm_sm2_get_base_affine_mont(&Gm);
    fixed_scalar(d, 0x12);

    for (int rr = 0; rr < REPEAT_MUL; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 20; i++) {
            asm_sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Gm, d);
            sink_sparse(&R, sizeof(R), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            T = R;
            sink_sparse(&T, sizeof(T), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            asm_sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Gm, d);
            sink_sparse(&R, sizeof(R), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_MUL
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_MUL
        );

        sink_mem(&R, sizeof(R));
    }

    return make_result(cyc, ns, REPEAT_MUL);
}
/* ========= gmssl point_mul_G ========= */

static bench_result_t bench_gmssl_point_mul_generator(void) {
    SM2_Z256_POINT R, T;
    sm2_z256_t d;
    uint8_t d_bytes[32];
    uint64_t cyc[REPEAT_MUL], ns[REPEAT_MUL];

    fixed_scalar(d_bytes, 0x12);
    z256_from_be32(d, d_bytes);

    for (int rr = 0; rr < REPEAT_MUL; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 20; i++) {
            sm2_z256_point_mul_generator(&R, d);
            sink_sparse(&R, sizeof(R), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            memcpy(&T, &R, sizeof(T));
            sink_sparse(&T, sizeof(T), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            sm2_z256_point_mul_generator(&R, d);
            sink_sparse(&R, sizeof(R), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_MUL
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_MUL
        );

        sink_mem(&R, sizeof(R));
    }

    return make_result(cyc, ns, REPEAT_MUL);
}

/* ========= self-c point_mul_P ========= */

static bench_result_t bench_self_c_point_mul(void) {
    sm2_affine_t Gm, Pm;
    sm2_jacobian_t R, T;
    uint8_t k[32], p_scalar[32];
    uint64_t cyc[REPEAT_MUL], ns[REPEAT_MUL];

    c_sm2_get_base_affine_mont(&Gm);

    fixed_scalar(p_scalar, 0x23);
    fixed_scalar(k, 0x45);

    /* 固定构造点 Pm = [p_scalar]G，保持 Mont affine */
    c_sm2_scalar_mul_window_ct_schemeB_mont_to_affine(&Pm, &Gm, p_scalar);

    for (int rr = 0; rr < REPEAT_MUL; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 20; i++) {
            c_sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Pm, k);
            sink_sparse(&R, sizeof(R), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            T = R;
            sink_sparse(&T, sizeof(T), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            c_sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Pm, k);
            sink_sparse(&R, sizeof(R), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_MUL
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_MUL
        );

        sink_mem(&R, sizeof(R));
    }

    return make_result(cyc, ns, REPEAT_MUL);
}

/* ========= self-asm point_mul_P ========= */

static bench_result_t bench_self_asm_point_mul(void) {
    sm2_affine_t Gm, Pm;
    sm2_jacobian_t R, T;
    uint8_t k[32], p_scalar[32];
    uint64_t cyc[REPEAT_MUL], ns[REPEAT_MUL];

    asm_sm2_get_base_affine_mont(&Gm);

    fixed_scalar(p_scalar, 0x23);
    fixed_scalar(k, 0x45);

    /* 固定构造点 Pm = [p_scalar]G，保持 Mont affine */
    asm_sm2_scalar_mul_window_ct_schemeB_mont_to_affine(&Pm, &Gm, p_scalar);

    for (int rr = 0; rr < REPEAT_MUL; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 20; i++) {
            asm_sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Pm, k);
            sink_sparse(&R, sizeof(R), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            T = R;
            sink_sparse(&T, sizeof(T), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            asm_sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Pm, k);
            sink_sparse(&R, sizeof(R), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_MUL
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_MUL
        );

        sink_mem(&R, sizeof(R));
    }

    return make_result(cyc, ns, REPEAT_MUL);
}

/* ========= gmssl point_mul_P ========= */

static bench_result_t bench_gmssl_point_mul(void) {
    SM2_Z256_POINT R, T, pub;
    sm2_z256_t k, d;
    uint8_t k_bytes[32], d_bytes[32];
    uint64_t cyc[REPEAT_MUL], ns[REPEAT_MUL];

    /* 与 self 版对齐：先固定标量，构造固定点 P = [d]G */
    fixed_scalar(d_bytes, 0x23);
    z256_from_be32(d, d_bytes);
    sm2_z256_point_mul_generator(&pub, d);

    fixed_scalar(k_bytes, 0x45);
    z256_from_be32(k, k_bytes);

    for (int rr = 0; rr < REPEAT_MUL; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 20; i++) {
            sm2_z256_point_mul(&R, k, &pub);
            sink_sparse(&R, sizeof(R), i);
        }

        uint64_t c0 = rdcycle();
        uint64_t t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            memcpy(&T, &R, sizeof(T));
            sink_sparse(&T, sizeof(T), i);
        }
        uint64_t c1 = rdcycle();
        uint64_t t1 = now_ns();
        empty_cyc = c1 - c0;
        empty_ns  = t1 - t0;

        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_MUL; i++) {
            sm2_z256_point_mul(&R, k, &pub);
            sink_sparse(&R, sizeof(R), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(
            full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc,
            ITER_MUL
        );
        ns[rr] = div_round_down(
            full_ns > empty_ns ? (full_ns - empty_ns) : full_ns,
            ITER_MUL
        );

        sink_mem(&R, sizeof(R));
    }

    return make_result(cyc, ns, REPEAT_MUL);
}


int main(void) {
    printf("Unified benchmark (c vs asm vs gmssl)\n");
    printf("=====================================\n");
    printf("FP: repeat=%d rounds=%d\n", REPEAT_FP, ITER_FP);
    printf("INV: repeat=%d rounds=%d\n", REPEAT_FP, ITER_INV);
    printf("FN: repeat=%d rounds=%d\n", REPEAT_FN, ITER_FN);

    bench_result_t c_fp_mul      = bench_self_c_fp_mont_mul();
    bench_result_t asm_fp_mul    = bench_self_asm_fp_mont_mul();
    bench_result_t gmssl_fp_mul  = bench_gmssl_modp_mont_mul();

    bench_result_t c_fp_sqr      = bench_self_c_fp_mont_sqr();
    bench_result_t asm_fp_sqr    = bench_self_asm_fp_mont_sqr();
    bench_result_t gmssl_fp_sqr  = bench_gmssl_modp_mont_sqr();

    bench_result_t c_fp_add      = bench_self_c_fp_mont_add();
    bench_result_t asm_fp_add    = bench_self_asm_fp_mont_add();
    bench_result_t gmssl_fp_add  = bench_gmssl_modp_add();
    bench_precise_result_t c_fp_add_precise     = bench_self_c_fp_mont_add_precise();
    bench_precise_result_t asm_fp_add_precise   = bench_self_asm_fp_mont_add_precise();
    bench_precise_result_t gmssl_fp_add_precise = bench_gmssl_modp_add_precise();

    bench_result_t c_fp_sub      = bench_self_c_fp_mont_sub();
    bench_result_t asm_fp_sub    = bench_self_asm_fp_mont_sub();
    bench_result_t gmssl_fp_sub  = bench_gmssl_modp_sub();

    bench_result_t c_fp_to       = bench_self_c_fp_to_mont();
    bench_result_t asm_fp_to     = bench_self_asm_fp_to_mont();
    bench_result_t gmssl_fp_to   = bench_gmssl_modp_to_mont();

    bench_result_t c_fp_from     = bench_self_c_fp_from_mont();
    bench_result_t asm_fp_from   = bench_self_asm_fp_from_mont();
    bench_result_t gmssl_fp_from = bench_gmssl_modp_from_mont();

    bench_result_t c_fp_inv      = bench_self_c_fp_mont_inv();
    bench_result_t asm_fp_inv    = bench_self_asm_fp_mont_inv();
    bench_result_t gmssl_fp_inv  = bench_gmssl_modp_mont_inv();

    bench_result_t c_fn_mont     = bench_self_c_fn_mont_mul();
    bench_result_t asm_fn_mont   = bench_self_asm_fn_mont_mul();
    bench_result_t gmssl_fn_mont = bench_gmssl_modn_mont_mul();

    bench_result_t c_fn_raw      = bench_self_c_fn_mul();
    bench_result_t asm_fn_raw    = bench_self_asm_fn_mul();
    bench_result_t gmssl_fn_raw  = bench_gmssl_modn_mul();

    bench_result_t c_pm_g      = bench_self_c_point_mul_generator();
    bench_result_t asm_pm_g    = bench_self_asm_point_mul_generator();
    bench_result_t gmssl_pm_g  = bench_gmssl_point_mul_generator();

    bench_result_t c_pm_p      = bench_self_c_point_mul();
    bench_result_t asm_pm_p    = bench_self_asm_point_mul();
    bench_result_t gmssl_pm_p  = bench_gmssl_point_mul();

    printf("\n[Fp]\n");
    print_row_3way("fp_mont_mul",  c_fp_mul,  asm_fp_mul,  gmssl_fp_mul);
    print_row_3way("fp_mont_sqr",  c_fp_sqr,  asm_fp_sqr,  gmssl_fp_sqr);
    print_row_3way("fp_mont_add",  c_fp_add,  asm_fp_add,  gmssl_fp_add);
    print_row_3way("fp_mont_sub",  c_fp_sub,  asm_fp_sub,  gmssl_fp_sub);
    print_row_3way("fp_to_mont",   c_fp_to,   asm_fp_to,   gmssl_fp_to);
    print_row_3way("fp_from_mont", c_fp_from, asm_fp_from, gmssl_fp_from);
    print_row_3way("fp_mont_inv",  c_fp_inv,  asm_fp_inv,  gmssl_fp_inv);

    printf("\n[Fn]\n");
    print_row_3way("fn_mont_mul",  c_fn_mont, asm_fn_mont, gmssl_fn_mont);
    print_row_3way("fn_mul",       c_fn_raw,  asm_fn_raw,  gmssl_fn_raw);

    printf("\n[ScalarMul]\n");
    print_row_3way("point_mul_G", c_pm_g, asm_pm_g, gmssl_pm_g);
    print_row_3way("point_mul_P", c_pm_p, asm_pm_p, gmssl_pm_p);

    printf("\n[Fp Precise]\n");
    print_row_3way_precise("fp_mont_add", c_fp_add_precise, asm_fp_add_precise, gmssl_fp_add_precise);

    printf("\n[sink] %u\n", bench_sink);
    return 0;
}
