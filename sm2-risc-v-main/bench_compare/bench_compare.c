#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>

/* ===== （纯 C，Montgomery 域）===== */
#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/fp.h"
#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/fn.h"

#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/sm2_curve.h"
#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/sm2_scalar.h"
/* ===== GmSSL ===== */
#include "/home/sm2/GmSSL/include/gmssl/sm2.h"
#include "/home/sm2/GmSSL/include/gmssl/sm2_z256.h"

#define SINK_STRIDE 128
#ifndef BENCH_USE_ASM
#define BENCH_USE_ASM 0
#endif

#if BENCH_USE_ASM
#define SELF_NAME "self-asm"
#else
#define SELF_NAME "self"
#endif


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

volatile uint8_t bench_sink = 0;
bench_result_t zero = {0};


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

static bench_result_t make_result(const uint64_t* cyc, const uint64_t* ns, int n) {
	bench_result_t r;

	uint64_t cyc_sum = 0, cyc_min = cyc[0], cyc_max = cyc[0];
	uint64_t ns_sum = 0, ns_min = ns[0], ns_max = ns[0];

	for (int i = 0; i < n; i++) {
		cyc_sum += cyc[i];
		ns_sum += ns[i];

		if (cyc[i] < cyc_min) cyc_min = cyc[i];
		if (cyc[i] > cyc_max) cyc_max = cyc[i];

		if (ns[i] < ns_min) ns_min = ns[i];
		if (ns[i] > ns_max) ns_max = ns[i];
	}

	r.cyc_mean = cyc_sum / (uint64_t)n;
	r.cyc_min = cyc_min;
	r.cyc_max = cyc_max;

	r.ns_mean = ns_sum / (uint64_t)n;
	r.ns_min = ns_min;
	r.ns_max = ns_max;

	return r;
}




static void print_section(const char* title) {
	printf("\n[%s]\n", title);
}

static void print_row_2way(const char* name, bench_result_t self_r, bench_result_t gmssl_r) {
	double speedup = (double)gmssl_r.cyc_mean / (double)self_r.cyc_mean;
	printf("%-16s  %-8s=%-9llu gmssl=%-9llu speedup=%.2fx\n",
		name,
		SELF_NAME,
		(unsigned long long)self_r.cyc_mean,
		(unsigned long long)gmssl_r.cyc_mean,
		speedup);
}

static void fixed_scalar(uint8_t k[32], uint8_t v) {
    memset(k, v, 32);
    k[31] |= 1;
}

static void fixed_scalar64(sm2_z256_t k, uint64_t v) {
    memset(k, 0, sizeof(sm2_z256_t));
    k[0] = v | 1ULL;
}

static inline uint64_t div_round_down(uint64_t x, uint64_t y) {
    return y ? (x / y) : 0;
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

/* ====================  fp_mont_mul ==================== */
static bench_result_t bench_self_fp_mont_mul(void) {
    fp_t a, b, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    fp_to_mont(&a, &a);
    fp_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        /* warmup */
        for (int i = 0; i < 200; i++) {
            fp_mont_mul(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        /* empty loop: 保留赋值链和 sink 开销 */
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

        /* full loop */
        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            fp_mont_mul(&r, &aa, &bb);
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

/* ==================== GmSSL modp_mont_mul ==================== */
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

        /* warmup */
        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_mont_mul(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        /* empty loop: 保留 memcpy + sink 开销 */
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

        /* full loop */
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


/* ==================== fp_mont_sqr ==================== */
static bench_result_t bench_self_fp_mont_sqr(void) {
    fp_t a, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;

    fp_to_mont(&a, &a);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        /* warmup */
        for (int i = 0; i < 200; i++) {
            fp_mont_sqr(&r, &aa);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        /* empty loop: 保留赋值链 + sparse sink 开销 */
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

        /* full loop */
        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            fp_mont_sqr(&r, &aa);
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

/* ==================== GmSSL modp_mont_sqr ==================== */
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

        /* warmup */
        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_mont_sqr(r, aa);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        /* empty loop: 保留 memcpy + sparse sink 开销 */
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

        /* full loop */
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

/* ==================== fp_mont_inv ==================== */
static bench_result_t bench_self_fp_mont_inv(void) {
    fp_t a, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;

    fp_to_mont(&a, &a);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        /* warmup: inv 比较重，预热次数不用太大 */
        for (int i = 0; i < 20; i++) {
            fp_mont_inv(&r, &aa);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        /* empty loop: 保留赋值链 + sparse sink 开销 */
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

        /* full loop */
        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_INV; i++) {
            fp_mont_inv(&r, &aa);
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

/* ==================== GmSSL modp_mont_inv ==================== */
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

        /* warmup */
        for (int i = 0; i < 20; i++) {
            sm2_z256_modp_mont_inv(r, aa);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        /* empty loop: 保留 memcpy + sparse sink 开销 */
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

        /* full loop */
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

/* ==================== 你的 fp_mont_add ==================== */
static bench_result_t bench_self_fp_mont_add(void) {
    fp_t a, b, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    fp_to_mont(&a, &a);
    fp_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        /* warmup */
        for (int i = 0; i < 200; i++) {
            fp_mont_add(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        /* empty loop: 保留赋值链 + sparse sink 开销 */
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

        /* full loop */
        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            fp_mont_add(&r, &aa, &bb);
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


/* ==================== GmSSL modp_add ==================== */
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

        /* warmup */
        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_add(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        /* empty loop: 保留 memcpy + sparse sink 开销 */
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

        /* full loop */
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

/* ====================  fp_mont_sub ==================== */
static bench_result_t bench_self_fp_mont_sub(void) {
    fp_t a, b, r;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x56, sizeof(a));
    memset(&b, 0x12, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    fp_to_mont(&a, &a);
    fp_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        fp_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        /* warmup */
        for (int i = 0; i < 200; i++) {
            fp_mont_sub(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        /* empty loop: 保留赋值链 + sparse sink 开销 */
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

        /* full loop */
        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            fp_mont_sub(&r, &aa, &bb);
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


/* ==================== GmSSL modp_sub ==================== */
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

        /* warmup */
        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_sub(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        /* empty loop: 保留 memcpy + sparse sink 开销 */
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

        /* full loop */
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


/* ====================  fp_to_mont ==================== */
static bench_result_t bench_self_fp_to_mont(void) {
    fp_t a, r, tmp;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        /* warmup: 固定普通域输入 */
        for (int i = 0; i < 200; i++) {
            fp_to_mont(&r, &a);
            sink_sparse(&r, sizeof(r), i);
        }

        /* empty loop: 保留读输入 + sparse sink 开销 */
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

        /* full loop: 始终测 a -> mont(r) */
        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            fp_to_mont(&r, &a);
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

/* ==================== GmSSL modp_to_mont ==================== */
static bench_result_t bench_gmssl_modp_to_mont(void) {
    sm2_z256_t a, r, tmp;
    uint8_t a_bytes[32];
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    fixed_scalar(a_bytes, 0x12);
    z256_from_be32(a, a_bytes);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        /* warmup: 固定普通域输入 */
        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_to_mont(a, r);
            sink_sparse(r, sizeof(r), i);
        }

        /* empty loop: 保留读输入 + sparse sink 开销 */
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

        /* full loop: 始终测 a -> mont(r) */
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

/* ====================  fp_from_mont ==================== */
static bench_result_t bench_self_fp_from_mont(void) {
    fp_t a, r, tmp;
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    memset(&a, 0x12, sizeof(a));
    a.v[0] |= 1;
    fp_to_mont(&a, &a);   /* 固定 Montgomery 域输入 */

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        /* warmup */
        for (int i = 0; i < 200; i++) {
            fp_from_mont(&r, &a);
            sink_sparse(&r, sizeof(r), i);
        }

        /* empty loop: 保留读输入 + sparse sink 开销 */
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

        /* full loop: 始终测 mont(a) -> normal(r) */
        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FP; i++) {
            fp_from_mont(&r, &a);
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


/* ==================== GmSSL modp_from_mont ==================== */
static bench_result_t bench_gmssl_modp_from_mont(void) {
    sm2_z256_t a, r, tmp;
    uint8_t a_bytes[32];
    uint64_t cyc[REPEAT_FP], ns[REPEAT_FP];

    fixed_scalar(a_bytes, 0x12);
    z256_from_be32(a, a_bytes);
    sm2_z256_modp_to_mont(a, a);

    for (int rr = 0; rr < REPEAT_FP; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        /* warmup */
        for (int i = 0; i < 200; i++) {
            sm2_z256_modp_from_mont(r, a);
            sink_sparse(r, sizeof(r), i);
        }

        /* empty loop: 保留读输入 + sparse sink 开销 */
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

        /* full loop: 始终测 mont(a) -> normal(r) */
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


/* ====================  fn_mont_mul_ct ==================== */
static bench_result_t bench_self_fn_mont_mul(void) {
    fn_t a, b, r;
    uint64_t cyc[REPEAT_FN], ns[REPEAT_FN];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    /* 转进 Montgomery 域 */
    fn_to_mont(&a, &a);
    fn_to_mont(&b, &b);

    for (int rr = 0; rr < REPEAT_FN; rr++) {
        fn_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        /* warmup */
        for (int i = 0; i < 200; i++) {
            fn_mont_mul_ct(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        /* empty loop: 保留赋值链 + sparse sink 开销 */
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

        /* full loop */
        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            fn_mont_mul_ct(&r, &aa, &bb);
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

/* ==================== GmSSL modn_mont_mul ==================== */
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

        /* warmup */
        for (int i = 0; i < 200; i++) {
            sm2_z256_modn_mont_mul(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        /* empty loop: 保留 memcpy + sparse sink 开销 */
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

        /* full loop */
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


/* ====================  fn_mul ==================== */
static bench_result_t bench_self_fn_mul(void) {
    fn_t a, b, r;
    uint64_t cyc[REPEAT_FN], ns[REPEAT_FN];

    memset(&a, 0x12, sizeof(a));
    memset(&b, 0x34, sizeof(b));
    a.v[0] |= 1;
    b.v[0] |= 1;

    for (int rr = 0; rr < REPEAT_FN; rr++) {
        fn_t aa = a, bb = b, tmp;
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        /* warmup */
        for (int i = 0; i < 200; i++) {
            fn_mul(&r, &aa, &bb);
            sink_sparse(&r, sizeof(r), i);
            aa = r;
        }

        /* empty loop: 保留赋值链 + sparse sink 开销 */
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

        /* full loop */
        c0 = rdcycle();
        t0 = now_ns();
        for (int i = 0; i < ITER_FN; i++) {
            fn_mul(&r, &aa, &bb);
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

/* ==================== GmSSL modn_mul ==================== */
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

        /* warmup */
        for (int i = 0; i < 200; i++) {
            sm2_z256_modn_mul(r, aa, bb);
            sink_sparse(r, sizeof(r), i);
            memcpy(aa, r, sizeof(r));
        }

        /* empty loop: 保留 memcpy + sparse sink 开销 */
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

        /* full loop */
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

/* ====================  point_mul_generator ([d]G), Mont ====================*/
static bench_result_t bench_self_point_mul_generator(void) {
    sm2_affine_t Gm;
    sm2_jacobian_t R, T;
    uint8_t d[32];
    uint64_t cyc[REPEAT_MUL], ns[REPEAT_MUL];

    sm2_get_base_affine_mont(&Gm);
    fixed_scalar(d, 0x12);

    for (int rr = 0; rr < REPEAT_MUL; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 20; i++) {
            sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Gm, d);
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
            sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Gm, d);
            sink_sparse(&R, sizeof(R), i);
        }
        c1 = rdcycle();
        t1 = now_ns();
        full_cyc = c1 - c0;
        full_ns  = t1 - t0;

        cyc[rr] = div_round_down(full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc, ITER_MUL);
        ns[rr]  = div_round_down(full_ns  > empty_ns  ? (full_ns  - empty_ns)  : full_ns,  ITER_MUL);

        sink_mem(&R, sizeof(R));
    }

    return make_result(cyc, ns, REPEAT_MUL);
}

/* ====================  gmssl_point_mul_generator ([d]G), Mont ====================*/
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

        cyc[rr] = div_round_down(full_cyc > empty_cyc ? (full_cyc - empty_cyc) : full_cyc, ITER_MUL);
        ns[rr]  = div_round_down(full_ns  > empty_ns  ? (full_ns  - empty_ns)  : full_ns,  ITER_MUL);

        sink_mem(&R, sizeof(R));
    }

    return make_result(cyc, ns, REPEAT_MUL);
}

/* ====================  point_mul ([d]P), Mont ==================== */
static bench_result_t bench_self_point_mul(void) {
    sm2_affine_t Gm, Pm;
    sm2_jacobian_t R, T;
    uint8_t k[32], p_scalar[32];
    uint64_t cyc[REPEAT_MUL], ns[REPEAT_MUL];

    sm2_get_base_affine_mont(&Gm);

    fixed_scalar(p_scalar, 0x23);
    fixed_scalar(k, 0x45);

    /* 固定构造点 Pm = [p_scalar]G，保持 Mont affine */
    sm2_scalar_mul_window_ct_schemeB_mont_to_affine(&Pm, &Gm, p_scalar);

    for (int rr = 0; rr < REPEAT_MUL; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 20; i++) {
            sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Pm, k);
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
            sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&R, &Pm, k);
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

/* ==================== GmSSL point_mul ([d]P) ==================== */
static bench_result_t bench_gmssl_point_mul(void) {
    SM2_Z256_POINT R, T, pub;
    sm2_z256_t k, d;
    uint8_t k_bytes[32], d_bytes[32];
    uint64_t cyc[REPEAT_MUL], ns[REPEAT_MUL];

    /* FIX: 与 self 版对齐，用固定标量 0x23 生成相同的 P */
    fixed_scalar(d_bytes, 0x23);
    z256_from_be32(d, d_bytes);
    sm2_z256_point_mul_generator(&pub, d);   /* pub = [0x23...23]G，和 self 的 Pm 等价 */

    fixed_scalar(k_bytes, 0x45);
    z256_from_be32(k, k_bytes);

    for (int rr = 0; rr < REPEAT_MUL; rr++) {
        uint64_t empty_cyc, empty_ns, full_cyc, full_ns;

        for (int i = 0; i < 20; i++) {
            sm2_z256_point_mul(&R, k, &pub);   /* FIX: 用固定 pub 替代随机 kb.public_key */
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
	printf("Compare benchmark (%s vs gmssl, cycles)\n", SELF_NAME);
	printf("==================================\n");
	printf("FP   : repeat=%d rounds=%d\n", REPEAT_FP, ITER_FP);
    printf("INV  : repeat=%d rounds=%d\n", REPEAT_FP, ITER_INV);
    printf("FN   : repeat=%d rounds=%d\n", REPEAT_FN, ITER_FN);
    printf("MUL  : repeat=%d rounds=%d\n", REPEAT_MUL, ITER_MUL);

	bench_result_t self_fp_mul = bench_self_fp_mont_mul();
	bench_result_t gmssl_fp_mul = bench_gmssl_modp_mont_mul();

	bench_result_t self_fp_sqr = bench_self_fp_mont_sqr();
	bench_result_t gmssl_fp_sqr = bench_gmssl_modp_mont_sqr();

	bench_result_t self_fp_inv = bench_self_fp_mont_inv();
	bench_result_t gmssl_fp_inv = bench_gmssl_modp_mont_inv();

	bench_result_t self_fp_add = bench_self_fp_mont_add();
	bench_result_t gmssl_fp_add = bench_gmssl_modp_add();

	bench_result_t self_fp_sub = bench_self_fp_mont_sub();
	bench_result_t gmssl_fp_sub = bench_gmssl_modp_sub();

	bench_result_t self_fp_to = bench_self_fp_to_mont();
	bench_result_t gmssl_fp_to = bench_gmssl_modp_to_mont();

	bench_result_t self_fp_from = bench_self_fp_from_mont();
	bench_result_t gmssl_fp_from = bench_gmssl_modp_from_mont();

	bench_result_t self_fn_raw = bench_self_fn_mont_mul();
	bench_result_t gmssl_fn_raw = bench_gmssl_modn_mont_mul();

	bench_result_t self_fn_full = bench_self_fn_mul();
	bench_result_t gmssl_fn_full = bench_gmssl_modn_mul();

	bench_result_t self_pm_g = bench_self_point_mul_generator();
	bench_result_t gmssl_pm_g = bench_gmssl_point_mul_generator();

	bench_result_t self_pm_p = bench_self_point_mul();
	bench_result_t gmssl_pm_p = bench_gmssl_point_mul();

	print_section("Fp");
	print_row_2way("fp_mont_mul", self_fp_mul, gmssl_fp_mul);
	print_row_2way("fp_mont_sqr", self_fp_sqr, gmssl_fp_sqr);
	print_row_2way("fp_mont_inv", self_fp_inv, gmssl_fp_inv);
	print_row_2way("fp_mont_add", self_fp_add, gmssl_fp_add);
	print_row_2way("fp_mont_sub", self_fp_sub, gmssl_fp_sub);
	print_row_2way("fp_to_mont", self_fp_to, gmssl_fp_to);
	print_row_2way("fp_from_mont", self_fp_from, gmssl_fp_from);

	print_section("Fn");
	print_row_2way("fn_mont_mul", self_fn_raw, gmssl_fn_raw);
	print_row_2way("fn_mul", self_fn_full, gmssl_fn_full);

	print_section("ScalarMul");
	print_row_2way("point_mul_G", self_pm_g, gmssl_pm_g);
	print_row_2way("point_mul_P", self_pm_p, gmssl_pm_p);

	printf("\n[sink] %u\n", bench_sink);
	return 0;
}