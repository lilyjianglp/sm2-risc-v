/*
 * Focused fp_mont_inv benchmark.
 *
 * This keeps the measurement surface small while tuning Pornin full ASM.
 */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "fp.h"

#ifdef PORNIN_FULL_INV_PROFILE
extern uint64_t pornin_full_prof_count;
extern uint64_t pornin_full_prof_extract;
extern uint64_t pornin_full_prof_inner31;
extern uint64_t pornin_full_prof_update_ab;
extern uint64_t pornin_full_prof_update_uv;
extern uint64_t pornin_full_prof_tail14;
extern uint64_t pornin_full_prof_tail_uv;
extern uint64_t pornin_full_prof_final_mul;
#endif

#ifndef INV_PERF_ITERS
#define INV_PERF_ITERS 10000
#endif

#ifndef INV_PERF_REPEAT
#define INV_PERF_REPEAT 9
#endif

#ifdef USE_PORNIN_FULL_INV_ASM
#define INV_LABEL "fp_mont_inv (Pornin full ASM)"
#elif defined(USE_PORNIN_INV) && defined(USE_PORNIN_INNER31_ASM) && defined(USE_PORNIN_UPDATE_UV_ASM)
#define INV_LABEL "fp_mont_inv (Pornin C+ASM)"
#elif defined(USE_PORNIN_INV)
#define INV_LABEL "fp_mont_inv (Pornin C)"
#else
#define INV_LABEL "fp_mont_inv (C)"
#endif

static volatile uint8_t bench_sink;

static inline uint64_t rdcycle(void)
{
    uint64_t c;
    asm volatile ("rdcycle %0" : "=r"(c));
    return c;
}

static inline void sink_mem(const void *p, size_t n)
{
    const volatile uint8_t *b = (const volatile uint8_t *)p;
    for (size_t i = 0; i < n; i++) {
        bench_sink ^= b[i];
    }
}

static void print_stats(const char *name, const uint64_t *v, int n)
{
    uint64_t sum = 0;
    uint64_t min = v[0];
    uint64_t max = v[0];

    for (int i = 0; i < n; i++) {
        sum += v[i];
        if (v[i] < min) {
            min = v[i];
        }
        if (v[i] > max) {
            max = v[i];
        }
    }

    printf("%-38s : mean=%" PRIu64 ", min=%" PRIu64 ", max=%" PRIu64 "\n",
           name, sum / (uint64_t)n, min, max);
}

#ifdef PORNIN_FULL_INV_PROFILE
static void print_profile_item(const char *name, uint64_t total, uint64_t count)
{
    printf("%-38s : total=%" PRIu64 ", per_inv=%" PRIu64 "\n",
           name, total, count ? total / count : 0);
}
#endif

int main(void)
{
    uint64_t cyc[INV_PERF_REPEAT];
    fp_t xi;
    fp_t ri;

    printf("Focused fp_mont_inv benchmark\n");
    printf("=============================\n");
    printf("iters=%d, repeat=%d\n\n", INV_PERF_ITERS, INV_PERF_REPEAT);

    memset(&xi, 0x9a, sizeof(xi));
    xi.v[0] |= 1;
    fp_to_mont(&xi, &xi);

    for (int rep = 0; rep < INV_PERF_REPEAT; rep++) {
        for (int i = 0; i < 10; i++) {
            fp_mont_inv(&ri, &xi);
            sink_mem(&ri, sizeof(ri));
        }

        uint64_t start = rdcycle();
        for (int i = 0; i < INV_PERF_ITERS; i++) {
            fp_mont_inv(&ri, &xi);
            fp_mont_add(&xi, &xi, &ri);
        }
        cyc[rep] = (rdcycle() - start) / (uint64_t)INV_PERF_ITERS;

        sink_mem(&ri, sizeof(ri));
        sink_mem(&xi, sizeof(xi));
    }

    print_stats(INV_LABEL, cyc, INV_PERF_REPEAT);

#ifdef PORNIN_FULL_INV_PROFILE
    printf("\nPornin full section profile\n");
    printf("---------------------------\n");
    printf("%-38s : %" PRIu64 "\n", "profiled invocations", pornin_full_prof_count);
    print_profile_item("extract x16", pornin_full_prof_extract, pornin_full_prof_count);
    print_profile_item("inner31 x16", pornin_full_prof_inner31, pornin_full_prof_count);
    print_profile_item("update_ab x16", pornin_full_prof_update_ab, pornin_full_prof_count);
    print_profile_item("update_uv x16", pornin_full_prof_update_uv, pornin_full_prof_count);
    print_profile_item("tail14", pornin_full_prof_tail14, pornin_full_prof_count);
    print_profile_item("tail_uv", pornin_full_prof_tail_uv, pornin_full_prof_count);
    print_profile_item("final fp_mont_mul", pornin_full_prof_final_mul, pornin_full_prof_count);
#endif

    printf("\n[sink] %u\n", (unsigned)bench_sink);
    return 0;
}
