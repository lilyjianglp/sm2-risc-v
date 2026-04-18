#include <stdio.h>
#include <string.h>
#include "fn.h"

static int test_fn_mul_correctness(void) {
    /* 已知输入输出的测试向量 */
    /* a = 1, b = 1 -> a*b mod n = 1 */
    fn_t a, b, r;
    fn_zero(&a); a.v[0] = 1;
    fn_zero(&b); b.v[0] = 1;
    fn_mul(&r, &a, &b);
    if (r.v[0] != 1 || r.v[1] || r.v[2] || r.v[3]) {
        printf("FAIL: 1*1 != 1\n");
        return 0;
    }

    /* a = 2, b = 3 -> 6 */
    fn_zero(&a); a.v[0] = 2;
    fn_zero(&b); b.v[0] = 3;
    fn_mul(&r, &a, &b);
    if (r.v[0] != 6 || r.v[1] || r.v[2] || r.v[3]) {
        printf("FAIL: 2*3 != 6, got %llx\n", (unsigned long long)r.v[0]);
        return 0;
    }

    /* a * (n-1) mod n = n - a，用 a=2 验证 */
    const fn_t *n = fn_modulus_n();
    fn_t n_minus_1 = *n;
    /* n-1: 直接减1 */
    n_minus_1.v[0] -= 1;
    fn_zero(&a); a.v[0] = 2;
    fn_mul(&r, &a, &n_minus_1);
    /* 期望结果 = n - 2 */
    fn_t expected = *n;
    expected.v[0] -= 2;
    if (memcmp(&r, &expected, sizeof(fn_t)) != 0) {
        printf("FAIL: 2*(n-1) != n-2\n");
        printf("  got:      %016llx %016llx %016llx %016llx\n",
               (unsigned long long)r.v[3], (unsigned long long)r.v[2],
               (unsigned long long)r.v[1], (unsigned long long)r.v[0]);
        printf("  expected: %016llx %016llx %016llx %016llx\n",
               (unsigned long long)expected.v[3], (unsigned long long)expected.v[2],
               (unsigned long long)expected.v[1], (unsigned long long)expected.v[0]);
        return 0;
    }

    /* 交换律：a*b == b*a */
    fn_zero(&a); a.v[0] = 0xDEADBEEFULL; a.v[1] = 0x12345678ULL;
    fn_zero(&b); b.v[0] = 0xCAFEBABEULL; b.v[1] = 0x87654321ULL;
    fn_t r2;
    fn_mul(&r,  &a, &b);
    fn_mul(&r2, &b, &a);
    if (memcmp(&r, &r2, sizeof(fn_t)) != 0) {
        printf("FAIL: a*b != b*a\n");
        return 0;
    }

    /* 结合律：(a*b)*c == a*(b*c) */
    fn_t c, r3, r4;
    fn_zero(&c); c.v[0] = 0x11223344ULL; c.v[2] = 0xAABBCCDDULL;
    fn_mul(&r3, &r, &c);   /* (a*b)*c */
    fn_mul(&r,  &b, &c);
    fn_mul(&r4, &a, &r);   /* a*(b*c) */
    if (memcmp(&r3, &r4, sizeof(fn_t)) != 0) {
        printf("FAIL: (a*b)*c != a*(b*c)\n");
        return 0;
    }

    printf("PASS: fn_mul correctness OK\n");
    return 1;
}

int main(void) {
    return test_fn_mul_correctness() ? 0 : 1;
}