#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "fp.h"

/* 待测汇编实现 */
extern void fp_mont_add(fp_t *r, const fp_t *a, const fp_t *b);

/* 来自 fp_montgomery.c */
extern const fp_t FP_P;

/* ============================================================
 * 工具函数
 * ============================================================ */

static inline uint64_t ct_mask_u64(uint64_t bit01) {
    return 0ULL - (bit01 & 1ULL);
}

static inline uint64_t addc_u64(uint64_t a, uint64_t b, uint64_t *carry) {
    uint64_t c  = *carry;
    uint64_t s  = a + b;
    uint64_t c1 = (s < a);
    uint64_t s2 = s + c;
    uint64_t c2 = (s2 < s);
    *carry = c1 | c2;
    return s2;
}

static inline uint64_t subb_u64(uint64_t a, uint64_t b, uint64_t *borrow) {
    uint64_t c  = *borrow;
    uint64_t d  = a - b;
    uint64_t b1 = (a < b);
    uint64_t d2 = d - c;
    uint64_t b2 = (d < c);
    *borrow = b1 | b2;
    return d2;
}

static int fp_equal(const fp_t *a, const fp_t *b) {
    return (a->v[0] == b->v[0] &&
            a->v[1] == b->v[1] &&
            a->v[2] == b->v[2] &&
            a->v[3] == b->v[3]);
}

static void fp_print(const char *name, const fp_t *a) {
    printf("%s = 0x%016llx_%016llx_%016llx_%016llx\n",
           name,
           (unsigned long long)a->v[3],
           (unsigned long long)a->v[2],
           (unsigned long long)a->v[1],
           (unsigned long long)a->v[0]);
}

static int fp_lt(const fp_t *a, const fp_t *b) {
    for (int i = 3; i >= 0; --i) {
        if (a->v[i] < b->v[i]) return 1;
        if (a->v[i] > b->v[i]) return 0;
    }
    return 0;
}

static int fp_lt_p(const fp_t *a) {
    return fp_lt(a, &FP_P);
}

static uint64_t rand64(void) {
    uint64_t x = 0;
    x ^= ((uint64_t)rand() & 0xFFFFu) <<  0;
    x ^= ((uint64_t)rand() & 0xFFFFu) << 16;
    x ^= ((uint64_t)rand() & 0xFFFFu) << 32;
    x ^= ((uint64_t)rand() & 0xFFFFu) << 48;
    return x;
}

static void fp_rand_mod_p(fp_t *x) {
    do {
        x->v[0] = rand64();
        x->v[1] = rand64();
        x->v[2] = rand64();
        x->v[3] = rand64();
    } while (!fp_lt_p(x));
}

static void fp_sub_small_from_p(fp_t *x, uint64_t k) {
    uint64_t borrow = 0;
    x->v[0] = subb_u64(FP_P.v[0], k, &borrow);
    x->v[1] = subb_u64(FP_P.v[1], 0, &borrow);
    x->v[2] = subb_u64(FP_P.v[2], 0, &borrow);
    x->v[3] = subb_u64(FP_P.v[3], 0, &borrow);
}

static void fp_from_u64(fp_t *x, uint64_t k) {
    x->v[0] = k;
    x->v[1] = 0;
    x->v[2] = 0;
    x->v[3] = 0;
}

/* ============================================================
 * C 参考实现
 * ============================================================ */

static void fp_mont_add_ref(fp_t *r, const fp_t *a, const fp_t *b) {
    uint64_t carry = 0;
    fp_t t;

    t.v[0] = addc_u64(a->v[0], b->v[0], &carry);
    t.v[1] = addc_u64(a->v[1], b->v[1], &carry);
    t.v[2] = addc_u64(a->v[2], b->v[2], &carry);
    t.v[3] = addc_u64(a->v[3], b->v[3], &carry);

    uint64_t borrow = 0;
    fp_t u;

    u.v[0] = subb_u64(t.v[0], FP_P.v[0], &borrow);
    u.v[1] = subb_u64(t.v[1], FP_P.v[1], &borrow);
    u.v[2] = subb_u64(t.v[2], FP_P.v[2], &borrow);
    u.v[3] = subb_u64(t.v[3], FP_P.v[3], &borrow);

    uint64_t m = ct_mask_u64((carry & 1ULL) | ((borrow ^ 1ULL) & 1ULL));

    r->v[0] = (u.v[0] & m) | (t.v[0] & ~m);
    r->v[1] = (u.v[1] & m) | (t.v[1] & ~m);
    r->v[2] = (u.v[2] & m) | (t.v[2] & ~m);
    r->v[3] = (u.v[3] & m) | (t.v[3] & ~m);
}

/* ============================================================
 * 对拍核心
 * ============================================================ */

static int g_case_id = 0;

static int report_fail(const char *tag,
                       const fp_t *a,
                       const fp_t *b,
                       const fp_t *ref,
                       const fp_t *got) {
    printf("FAIL at case #%d [%s]\n", g_case_id, tag);
    fp_print("a   ", a);
    fp_print("b   ", b);
    fp_print("ref ", ref);
    fp_print("asm ", got);
    return 0;
}

static int one_case_basic(const char *tag, const fp_t *a, const fp_t *b) {
    fp_t ref, got;
    fp_mont_add_ref(&ref, a, b);
    fp_mont_add(&got, a, b);
    g_case_id++;
    if (!fp_equal(&ref, &got)) {
        return report_fail(tag, a, b, &ref, &got);
    }
    return 1;
}

static int one_case_alias_a(const char *tag, const fp_t *a, const fp_t *b) {
    fp_t ref, got, ain;
    fp_mont_add_ref(&ref, a, b);
    ain = *a;
    fp_mont_add(&ain, &ain, b);   /* r == a */
    got = ain;
    g_case_id++;
    if (!fp_equal(&ref, &got)) {
        return report_fail(tag, a, b, &ref, &got);
    }
    return 1;
}

static int one_case_alias_b(const char *tag, const fp_t *a, const fp_t *b) {
    fp_t ref, got, bin;
    fp_mont_add_ref(&ref, a, b);
    bin = *b;
    fp_mont_add(&bin, a, &bin);   /* r == b */
    got = bin;
    g_case_id++;
    if (!fp_equal(&ref, &got)) {
        return report_fail(tag, a, b, &ref, &got);
    }
    return 1;
}

static int one_case_all(const char *tag, const fp_t *a, const fp_t *b) {
    if (!one_case_basic(tag, a, b)) return 0;
    if (!one_case_alias_a(tag, a, b)) return 0;
    if (!one_case_alias_b(tag, a, b)) return 0;
    return 1;
}

/* ============================================================
 * 主程序
 * ============================================================ */

int main(void) {
    srand(1);

    fp_t a, b;
    fp_t zero = {{0,0,0,0}};
    fp_t one  = {{1,0,0,0}};
    fp_t pm1  = {{
        0xFFFFFFFFFFFFFFFEULL,
        0xFFFFFFFF00000000ULL,
        0xFFFFFFFFFFFFFFFFULL,
        0xFFFFFFFEFFFFFFFFULL
    }};

    /* --------------------------------------------------------
     * 1) 基本边界测试
     * -------------------------------------------------------- */
    if (!one_case_all("basic: zero+zero", &zero, &zero)) return 1;
    if (!one_case_all("basic: zero+one",  &zero, &one )) return 1;
    if (!one_case_all("basic: one+one",   &one,  &one )) return 1;
    if (!one_case_all("basic: pm1+zero",  &pm1,  &zero)) return 1;
    if (!one_case_all("basic: pm1+one",   &pm1,  &one )) return 1;
    if (!one_case_all("basic: pm1+pm1",   &pm1,  &pm1 )) return 1;

    /* --------------------------------------------------------
     * 2) 单 limb / 连锁进位测试
     * -------------------------------------------------------- */
    a = (fp_t){{0xFFFFFFFFFFFFFFFFULL, 0, 0, 0}};
    b = (fp_t){{1,0,0,0}};
    if (!one_case_all("carry: limb0", &a, &b)) return 1;

    a = (fp_t){{0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0, 0}};
    b = (fp_t){{1,0,0,0}};
    if (!one_case_all("carry: limb0->1", &a, &b)) return 1;

    a = (fp_t){{0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
                0xFFFFFFFFFFFFFFFFULL, 0}};
    b = (fp_t){{1,0,0,0}};
    if (!one_case_all("carry: limb0->2", &a, &b)) return 1;

    a = (fp_t){{0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
                0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL}};
    b = (fp_t){{1,0,0,0}};
    if (!one_case_all("carry: full chain", &a, &b)) return 1;

    /* --------------------------------------------------------
     * 3) 近 p 的压力测试
     *    a = p-i, b = p-j
     * -------------------------------------------------------- */
    for (uint64_t i = 1; i <= 20000; ++i) {
        uint64_t j = i + 1;
        fp_sub_small_from_p(&a, i);
        fp_sub_small_from_p(&b, j);
        if (!one_case_all("near-p: p-i + p-j", &a, &b)) return 1;
    }

    /* --------------------------------------------------------
     * 4) 小数与近 p 混合
     * -------------------------------------------------------- */
    for (uint64_t i = 0; i <= 20000; ++i) {
        fp_from_u64(&a, i);
        fp_sub_small_from_p(&b, i + 1);
        if (!one_case_all("mixed: small + near-p", &a, &b)) return 1;
    }

    /* --------------------------------------------------------
     * 5) 只设某一 limb 的定向测试
     * -------------------------------------------------------- */
    for (uint64_t k = 0; k < 10000; ++k) {
        a = (fp_t){{rand64(), 0, 0, 0}};
        b = (fp_t){{rand64(), 0, 0, 0}};
        if (fp_lt_p(&a) && fp_lt_p(&b)) {
            if (!one_case_all("pattern: limb0 only", &a, &b)) return 1;
        }

        a = (fp_t){{0, rand64(), 0, 0}};
        b = (fp_t){{0, rand64(), 0, 0}};
        if (fp_lt_p(&a) && fp_lt_p(&b)) {
            if (!one_case_all("pattern: limb1 only", &a, &b)) return 1;
        }

        a = (fp_t){{0, 0, rand64(), 0}};
        b = (fp_t){{0, 0, rand64(), 0}};
        if (fp_lt_p(&a) && fp_lt_p(&b)) {
            if (!one_case_all("pattern: limb2 only", &a, &b)) return 1;
        }

        a = (fp_t){{0, 0, 0, rand64()}};
        b = (fp_t){{0, 0, 0, rand64()}};
        if (fp_lt_p(&a) && fp_lt_p(&b)) {
            if (!one_case_all("pattern: limb3 only", &a, &b)) return 1;
        }
    }

    /* --------------------------------------------------------
     * 6) 大量随机测试
     * -------------------------------------------------------- */
    for (int i = 0; i < 500000; ++i) {
        fp_rand_mod_p(&a);
        fp_rand_mod_p(&b);
        if (!one_case_all("random", &a, &b)) return 1;
    }

    printf("All strengthened fp_mont_add differential tests passed: %d cases\n",
           g_case_id);
    return 0;
}
