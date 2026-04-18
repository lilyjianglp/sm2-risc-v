/*
 * correctness_compare.c
 *
 * Cross-check this project's SM2 field, scalar, and point arithmetic
 * against GmSSL.
 *
 * Coverage:
 * - fp_to_mont / fp_from_mont
 * - fp_mont_add / sub / mul / sqr / inv
 * - fn_mont_mul / fn_mul
 * - point_mul_generator
 * - point_mul
 *
 * The comparison strategy is:
 * - Compare internal Montgomery-domain values directly when layouts match
 * - Compare canonical 32-byte outputs for scalar and point results
 * - Stop on the first mismatch and print enough data to debug it
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

/* ===== Local implementation ===== */
#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/fp.h"
#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/fn.h"
#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/sm2_curve.h"
#include "/home/sm2/sm2-risc-v-main/sm2-risc-v-main/src/sm2_scalar.h"

/* ===== GmSSL ===== */
#include "/home/sm2/GmSSL/include/gmssl/sm2.h"
#include "/home/sm2/GmSSL/include/gmssl/sm2_z256.h"

#define TEST_ROUNDS 1000
#define STRESS_ROUNDS 256

static int g_total = 0;
static int g_pass  = 0;
static int g_fail  = 0;

static const uint8_t SM2_N_BYTES[32] = {
    0xFF,0xFF,0xFF,0xFE,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0x72,0x03,0xDF,0x6B,0x21,0xC6,0x05,0x2B,
    0x53,0xBB,0xF4,0x09,0x39,0xD5,0x41,0x23
};

/* =========================================================
 * Utility helpers
 * ========================================================= */

static void print_hex_bytes(const char *name, const uint8_t *p, size_t n) {
    printf("%s = ", name);
    for (size_t i = 0; i < n; i++) {
        printf("%02x", p[i]);
    }
    printf("\n");
}

static void print_fp_hex(const char *name, const fp_t *a) {
    print_hex_bytes(name, (const uint8_t *)a, sizeof(fp_t));
}

static void print_fn_hex(const char *name, const fn_t *a) {
    print_hex_bytes(name, (const uint8_t *)a, sizeof(fn_t));
}

static void print_z256_hex(const char *name, const sm2_z256_t a) {
    print_hex_bytes(name, (const uint8_t *)a, sizeof(sm2_z256_t));
}

static void mark_result(const char *name, int ok) {
    g_total++;
    if (ok) {
        g_pass++;
        printf("[PASS] %s\n", name);
    } else {
        g_fail++;
        printf("[FAIL] %s\n", name);
    }
}

static void fill_pattern_u8(uint8_t *buf, size_t n, uint8_t seed) {
    for (size_t i = 0; i < n; i++) {
        buf[i] = (uint8_t)(seed + (uint8_t)(13 * i));
    }
}

static void make_test_fp(fp_t *a, uint8_t seed) {
    memset(a, 0, sizeof(*a));
    fill_pattern_u8((uint8_t *)a, sizeof(*a), seed);
    a->v[0] |= 1ULL;
}

static void make_test_fn(fn_t *a, uint8_t seed) {
    uint8_t be[32];
    fill_pattern_u8(be, sizeof(be), seed);
    be[31] |= 1;
    fn_from_be(a, be);
}


/* GmSSL sm2_z256_t and local fp_t/fn_t share the same 4x64-bit little-endian limb layout. */
static void fp_to_gmssl(sm2_z256_t out, const fp_t *in) {
    memcpy(out, in, sizeof(sm2_z256_t));
}

static void gmssl_to_fp(fp_t *out, const sm2_z256_t in) {
    memcpy(out, in, sizeof(sm2_z256_t));
}

static void fn_to_gmssl(sm2_z256_t out, const fn_t *in) {
    memcpy(out, in, sizeof(sm2_z256_t));
}

static void gmssl_to_fn(fn_t *out, const sm2_z256_t in) {
    memcpy(out, in, sizeof(sm2_z256_t));
}

static int fp_equal(const fp_t *a, const fp_t *b) {
    return memcmp(a, b, sizeof(fp_t)) == 0;
}

static uint64_t subb_u64_test(uint64_t a, uint64_t b, uint64_t *borrow) {
    uint64_t bb = b + *borrow;
    uint64_t bcarry = (bb < b);
    uint64_t r = a - bb;
    *borrow = (uint64_t)((a < bb) | bcarry);
    return r;
}

static void fp_normalize_once_for_compare(fp_t *r, const fp_t *a) {
    uint64_t borrow = 0;
    fp_t u;

    u.v[0] = subb_u64_test(a->v[0], FP_P.v[0], &borrow);
    u.v[1] = subb_u64_test(a->v[1], FP_P.v[1], &borrow);
    u.v[2] = subb_u64_test(a->v[2], FP_P.v[2], &borrow);
    u.v[3] = subb_u64_test(a->v[3], FP_P.v[3], &borrow);

    uint64_t keep_a = 0 - borrow;
    uint64_t keep_u = ~keep_a;

    r->v[0] = (a->v[0] & keep_a) | (u.v[0] & keep_u);
    r->v[1] = (a->v[1] & keep_a) | (u.v[1] & keep_u);
    r->v[2] = (a->v[2] & keep_a) | (u.v[2] & keep_u);
    r->v[3] = (a->v[3] & keep_a) | (u.v[3] & keep_u);
}

static int bytes_equal(const uint8_t *a, const uint8_t *b, size_t n) {
    return memcmp(a, b, n) == 0;
}

static int scalar_is_zero_bytes(const uint8_t x[32]) {
    uint8_t acc = 0;
    for (int i = 0; i < 32; i++) acc |= x[i];
    return acc == 0;
}

static int scalar_is_order_bytes(const uint8_t x[32]) {
    return bytes_equal(x, SM2_N_BYTES, 32);
}

static void print_bytes32(const char *name, const uint8_t x[32]) {
    printf("%s = ", name);
    for (int i = 0; i < 32; i++) {
        printf("%02x", x[i]);
    }
    printf("\n");
}

static void scalar_set_u64(uint8_t out[32], uint64_t x) {
    memset(out, 0, 32);
    for (int i = 31; i >= 24; i--) {
        out[i] = (uint8_t)x;
        x >>= 8;
    }
}

static void scalar_copy(uint8_t out[32], const uint8_t in[32]) {
    memcpy(out, in, 32);
}

static void scalar_add_small(uint8_t out[32], const uint8_t in[32], uint32_t add) {
    uint32_t carry = add;
    memcpy(out, in, 32);
    for (int i = 31; i >= 0 && carry != 0; i--) {
        uint32_t v = (uint32_t)out[i] + (carry & 0xffU);
        out[i] = (uint8_t)v;
        carry = (carry >> 8) + (v >> 8);
    }
}

static void scalar_sub_small(uint8_t out[32], const uint8_t in[32], uint32_t sub) {
    uint32_t borrow = sub;
    memcpy(out, in, 32);
    for (int i = 31; i >= 0 && borrow != 0; i--) {
        uint32_t s = borrow & 0xffU;
        uint32_t v = out[i];
        out[i] = (uint8_t)(v - s);
        borrow = (borrow >> 8) + (v < s);
    }
}

static void scalar_set_bit(uint8_t out[32], int bit) {
    memset(out, 0, 32);
    if (bit >= 0 && bit < 256) {
        out[31 - (bit >> 3)] = (uint8_t)(1U << (bit & 7));
    }
}

static void scalar_fill_pattern(uint8_t out[32], int pattern_id) {
    memset(out, 0, 32);
    switch (pattern_id) {
    case 0:
        memset(out, 0x55, 32);
        break;
    case 1:
        memset(out, 0xaa, 32);
        break;
    case 2:
        memset(out, 0xff, 32);
        break;
    case 3:
        for (int i = 0; i < 32; i++) out[i] = (uint8_t)(0x11U * (uint8_t)(i + 1));
        break;
    case 4:
        for (int i = 0; i < 32; i++) out[i] = (uint8_t)(0xf0U ^ (uint8_t)(i * 17));
        break;
    default:
        fill_pattern_u8(out, 32, (uint8_t)(0x31 + pattern_id));
        break;
    }
}

static void self_fn_mont_to_bytes32(uint8_t out[32], const fn_t *in_mont) {
    fn_t t;
    fn_from_mont_pub(&t, in_mont);
    fn_to_be(out, &t);
}

static void self_fn_to_bytes32(uint8_t out[32], const fn_t *in) {
    fn_to_be(out, in);
}

static int self_jacobian_mont_to_bytes32(
    uint8_t x[32], uint8_t y[32], const sm2_jacobian_t *P)
{
    sm2_affine_t a_mont, a_plain;

    if (sm2_jacobian_is_infinity(P)) {
        return 0;
    }

    sm2_jacobian_mont_to_affine_mont(&a_mont, P);
    sm2_affine_from_mont(&a_plain, &a_mont);

    if (a_plain.infinity) {
        return 0;
    }

    fp_to_bytes(x, &a_plain.x);
    fp_to_bytes(y, &a_plain.y);
    return 1;
}

static int gmssl_point_to_bytes32(
    uint8_t x[32], uint8_t y[32], const SM2_Z256_POINT *P)
{
    sm2_z256_t gx, gy;

    if (sm2_z256_point_get_xy(P, gx, gy) != 1) {
        return 0;
    }

    sm2_z256_to_bytes(gx, x);
    sm2_z256_to_bytes(gy, y);
    return 1;
}

/* =========================================================
 * 1. fp_to_mont / fp_from_mont
 * ========================================================= */

static int test_fp_to_from_mont_once(uint8_t seed) {
    fp_t a, self_to, self_back;
    sm2_z256_t ga, g_to, g_back;
    uint8_t a_be[32], self_back_be[32], gm_back_be[32];

    make_test_fp(&a, seed);

    fp_to_mont(&self_to, &a);
    fp_from_mont(&self_back, &self_to);

    fp_to_gmssl(ga, &a);
    sm2_z256_modp_to_mont(ga, g_to);
    sm2_z256_modp_from_mont(g_back, g_to);

    fp_to_bytes(a_be, &a);
    fp_to_bytes(self_back_be, &self_back);

    {
        fp_t tmp;
        gmssl_to_fp(&tmp, g_back);
        fp_to_bytes(gm_back_be, &tmp);
    }

    if (!bytes_equal(a_be, self_back_be, 32)) {
        printf("Self roundtrip failed, seed=%u\n", seed);
        print_bytes32("input", a_be);
        print_bytes32("self_back", self_back_be);
        return 0;
    }

    if (!bytes_equal(a_be, gm_back_be, 32)) {
        printf("GmSSL roundtrip failed, seed=%u\n", seed);
        print_bytes32("input", a_be);
        print_bytes32("gm_back", gm_back_be);
        return 0;
    }

    return 1;
}

static int compare_self_gmssl_points(
    const char *label,
    const sm2_jacobian_t *self_R,
    const SM2_Z256_POINT *gm_R,
    const uint8_t *d1,
    const uint8_t *d2)
{
    uint8_t self_x[32], self_y[32], gm_x[32], gm_y[32];
    int self_ok = self_jacobian_mont_to_bytes32(self_x, self_y, self_R);
    int gm_ok = gmssl_point_to_bytes32(gm_x, gm_y, gm_R);

    if (!self_ok && !gm_ok) {
        return 1;
    }

    if (self_ok != gm_ok ||
        !bytes_equal(self_x, gm_x, 32) ||
        !bytes_equal(self_y, gm_y, 32)) {
        printf("Mismatch in %s\n", label);
        if (d1) print_bytes32("scalar_1", d1);
        if (d2) print_bytes32("scalar_2", d2);
        printf("self_export_ok = %d\n", self_ok);
        printf("gmssl_export_ok = %d\n", gm_ok);
        if (self_ok) {
            print_bytes32("self_x", self_x);
            print_bytes32("self_y", self_y);
        }
        if (gm_ok) {
            print_bytes32("gm_x", gm_x);
            print_bytes32("gm_y", gm_y);
        }
        return 0;
    }

    return 1;
}

static void test_fp_to_from_mont(void) {
    int ok = 1;
    for (int i = 0; i < TEST_ROUNDS; i++) {
        if (!test_fp_to_from_mont_once((uint8_t)(i + 1))) {
            ok = 0;
            break;
        }
    }
    mark_result("fp_to_mont + fp_from_mont", ok);
}

/* =========================================================
 * 2. fp add/sub/mul/sqr/inv
 * ========================================================= */

static int test_fp_add_once(uint8_t s1, uint8_t s2) {
    fp_t a, b, a_m, b_m, self_r;
    sm2_z256_t ga, gb, gr;

    make_test_fp(&a, s1);
    make_test_fp(&b, s2);

    fp_to_mont(&a_m, &a);
    fp_to_mont(&b_m, &b);
    fp_mont_add(&self_r, &a_m, &b_m);

    fp_to_gmssl(ga, &a_m);
    fp_to_gmssl(gb, &b_m);
    sm2_z256_modp_add(gr, ga, gb);

    {
        fp_t gm_r;
        gmssl_to_fp(&gm_r, gr);

        fp_t self_norm;
        fp_normalize_once_for_compare(&self_norm, &self_r);

        if (!fp_equal(&self_norm, &gm_r)) {
            printf("Mismatch in fp_mont_add, seeds=(%u,%u)\n", s1, s2);
            print_fp_hex("self_r", &self_r);
            print_fp_hex("self_norm", &self_norm);
            print_fp_hex("gmssl_r", &gm_r);
            return 0;
        }
    }

    return 1;
}

static int test_fp_sub_once(uint8_t s1, uint8_t s2) {
    fp_t a, b, a_m, b_m, self_r;
    sm2_z256_t ga, gb, gr;

    make_test_fp(&a, s1);
    make_test_fp(&b, s2);

    fp_to_mont(&a_m, &a);
    fp_to_mont(&b_m, &b);
    fp_mont_sub(&self_r, &a_m, &b_m);

    fp_to_gmssl(ga, &a_m);
    fp_to_gmssl(gb, &b_m);
    sm2_z256_modp_sub(gr, ga, gb);

    {
        fp_t gm_r;
        gmssl_to_fp(&gm_r, gr);

        if (!fp_equal(&self_r, &gm_r)) {
            printf("Mismatch in fp_mont_sub, seeds=(%u,%u)\n", s1, s2);
            print_fp_hex("self_r", &self_r);
            print_fp_hex("gmssl_r", &gm_r);
            return 0;
        }
    }

    return 1;
}

static int test_fp_mul_once(uint8_t s1, uint8_t s2) {
    fp_t a, b, a_m, b_m, self_r;
    sm2_z256_t ga, gb, gr;

    make_test_fp(&a, s1);
    make_test_fp(&b, s2);

    fp_to_mont(&a_m, &a);
    fp_to_mont(&b_m, &b);
    fp_mont_mul(&self_r, &a_m, &b_m);

    fp_to_gmssl(ga, &a_m);
    fp_to_gmssl(gb, &b_m);
    sm2_z256_modp_mont_mul(gr, ga, gb);

    {
        fp_t gm_r;
        gmssl_to_fp(&gm_r, gr);

        if (!fp_equal(&self_r, &gm_r)) {
            printf("Mismatch in fp_mont_mul, seeds=(%u,%u)\n", s1, s2);
            print_fp_hex("self_r", &self_r);
            print_fp_hex("gmssl_r", &gm_r);
            return 0;
        }
    }

    return 1;
}

static int test_fp_sqr_once(uint8_t seed) {
    fp_t a, a_m, self_r, self_mul;
    sm2_z256_t ga, gr;

    make_test_fp(&a, seed);
    fp_to_mont(&a_m, &a);

    fp_mont_sqr(&self_r, &a_m);
    fp_mont_mul(&self_mul, &a_m, &a_m);

    if (!fp_equal(&self_r, &self_mul)) {
        printf("Internal mismatch: fp_mont_sqr != fp_mont_mul(a, a), seed=%u\n", seed);
        print_fp_hex("self_sqr", &self_r);
        print_fp_hex("self_mul", &self_mul);
        return 0;
    }

    fp_to_gmssl(ga, &a_m);
    sm2_z256_modp_mont_sqr(gr, ga);

    {
        fp_t gm_r;
        gmssl_to_fp(&gm_r, gr);

        if (!fp_equal(&self_r, &gm_r)) {
            printf("Mismatch in fp_mont_sqr, seed=%u\n", seed);
            print_fp_hex("self_r", &self_r);
            print_fp_hex("gmssl_r", &gm_r);
            return 0;
        }
    }

    return 1;
}

static int test_fp_inv_once(uint8_t seed) {
    fp_t a, a_m, self_inv, self_one;
    sm2_z256_t ga, gr;

    make_test_fp(&a, seed);
    fp_to_mont(&a_m, &a);

    fp_mont_inv(&self_inv, &a_m);
    fp_mont_mul(&self_one, &a_m, &self_inv);

    fp_to_gmssl(ga, &a_m);
    sm2_z256_modp_mont_inv(gr, ga);

    {
        fp_t gm_inv;
        gmssl_to_fp(&gm_inv, gr);

        if (!fp_equal(&self_inv, &gm_inv)) {
            printf("Mismatch in fp_mont_inv, seed=%u\n", seed);
            print_fp_hex("self_inv", &self_inv);
            print_fp_hex("gmssl_inv", &gm_inv);
            return 0;
        }
    }

    {
        fp_t mont_one_from_self, one_plain;
        memset(&one_plain, 0, sizeof(one_plain));
        one_plain.v[0] = 1;
        fp_to_mont(&mont_one_from_self, &one_plain);

        if (!fp_equal(&self_one, &mont_one_from_self)) {
            printf("Property check failed: a * inv(a) != 1 in Montgomery domain, seed=%u\n", seed);
            print_fp_hex("a*inv(a)", &self_one);
            print_fp_hex("mont_one", &mont_one_from_self);
            return 0;
        }
    }

    return 1;
}

static void test_fp_ops(void) {
    int ok_add = 1, ok_sub = 1, ok_mul = 1, ok_sqr = 1, ok_inv = 1;

    for (int i = 0; i < TEST_ROUNDS; i++) {
        uint8_t s1 = (uint8_t)(i + 1);
        uint8_t s2 = (uint8_t)(0x80u + i);

        if (ok_add && !test_fp_add_once(s1, s2)) ok_add = 0;
        if (ok_sub && !test_fp_sub_once(s1, s2)) ok_sub = 0;
        if (ok_mul && !test_fp_mul_once(s1, s2)) ok_mul = 0;
        if (ok_sqr && !test_fp_sqr_once(s1)) ok_sqr = 0;
        if (ok_inv && !test_fp_inv_once((uint8_t)(s1 | 1u))) ok_inv = 0;

        if (!ok_add || !ok_sub || !ok_mul || !ok_sqr || !ok_inv) {
            break;
        }
    }

    mark_result("fp_mont_add", ok_add);
    mark_result("fp_mont_sub", ok_sub);
    mark_result("fp_mont_mul", ok_mul);
    mark_result("fp_mont_sqr", ok_sqr);
    mark_result("fp_mont_inv", ok_inv);
}

/* =========================================================
 * 3. fn_mont_mul / fn_mul
 * ========================================================= */

static int test_fn_mont_mul_once(uint8_t s1, uint8_t s2) {
    fn_t a, b, a_m, b_m, self_r;
    sm2_z256_t ga, gb, gr;
    uint8_t self_be[32], gm_be[32];

    make_test_fn(&a, s1);
    make_test_fn(&b, s2);

    fn_to_mont_pub(&a_m, &a);
    fn_to_mont_pub(&b_m, &b);
    fn_mont_mul_ct(&self_r, &a_m, &b_m);
    self_fn_mont_to_bytes32(self_be, &self_r);

    fn_to_gmssl(ga, &a_m);
    fn_to_gmssl(gb, &b_m);
    sm2_z256_modn_mont_mul(gr, ga, gb);
    sm2_z256_modn_from_mont(gr, gr);
    sm2_z256_to_bytes(gr, gm_be);

    if (!bytes_equal(self_be, gm_be, 32)) {
        printf("Mismatch in fn_mont_mul, seeds=(%u,%u)\n", s1, s2);
        print_bytes32("self", self_be);
        print_bytes32("gmssl", gm_be);
        return 0;
    }

    return 1;
}

static int test_fn_mul_once(uint8_t s1, uint8_t s2) {
    fn_t a, b, self_r;
    sm2_z256_t ga, gb, gr;
    uint8_t self_be[32], gm_be[32];

    make_test_fn(&a, s1);
    make_test_fn(&b, s2);

    fn_mul(&self_r, &a, &b);
    self_fn_to_bytes32(self_be, &self_r);

    fn_to_gmssl(ga, &a);
    fn_to_gmssl(gb, &b);
    sm2_z256_modn_mul(gr, ga, gb);
    sm2_z256_to_bytes(gr, gm_be);

    if (!bytes_equal(self_be, gm_be, 32)) {
        printf("Mismatch in fn_mul, seeds=(%u,%u)\n", s1, s2);
        print_bytes32("self", self_be);
        print_bytes32("gmssl", gm_be);
        return 0;
    }

    return 1;
}

static void test_fn_ops(void) {
    int ok_mont = 1, ok_full = 1;

    for (int i = 0; i < TEST_ROUNDS; i++) {
        uint8_t s1 = (uint8_t)(i + 3);
        uint8_t s2 = (uint8_t)(0x55u + i);

        if (ok_mont && !test_fn_mont_mul_once(s1, s2)) ok_mont = 0;
        if (ok_full && !test_fn_mul_once(s1, s2)) ok_full = 0;

        if (!ok_mont || !ok_full) break;
    }

    mark_result("fn_mont_mul", ok_mont);
    mark_result("fn_mul", ok_full);
}

/* =========================================================
 * 4. point_mul_generator / point_mul
 * ========================================================= */

static int test_point_mul_generator_once(uint8_t seed) {
    sm2_affine_t Gm;
    sm2_jacobian_t self_R;
    uint8_t d[32];
    uint8_t self_x[32], self_y[32], gm_x[32], gm_y[32];

    SM2_Z256_POINT gm_R;
    sm2_z256_t gd;

    sm2_get_base_affine_mont(&Gm);

    memset(d, 0, sizeof(d));
    fill_pattern_u8(d, sizeof(d), seed);
    d[31] |= 1;

    sm2_z256_from_bytes(gd, d);

    sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&self_R, &Gm, d);
    sm2_z256_point_mul_generator(&gm_R, gd);

    if (!self_jacobian_mont_to_bytes32(self_x, self_y, &self_R)) {
        printf("self point_mul_generator export failed, seed=%u\n", seed);
        return 0;
    }

    if (!gmssl_point_to_bytes32(gm_x, gm_y, &gm_R)) {
        printf("gmssl point_mul_generator export failed, seed=%u\n", seed);
        return 0;
    }

    if (!bytes_equal(self_x, gm_x, 32) || !bytes_equal(self_y, gm_y, 32)) {
        printf("Mismatch in point_mul_generator, seed=%u\n", seed);
        print_bytes32("self_x", self_x);
        print_bytes32("self_y", self_y);
        print_bytes32("gm_x", gm_x);
        print_bytes32("gm_y", gm_y);
        return 0;
    }

    return 1;
}

static void test_point_mul_generator(void) {
    int ok = 1;
    for (int i = 0; i < TEST_ROUNDS; i++) {
        if (!test_point_mul_generator_once((uint8_t)(i + 7))) {
            ok = 0;
            break;
        }
    }
    mark_result("point_mul_generator", ok);
}


static int test_point_mul_once(uint8_t seed1, uint8_t seed2) {
    sm2_affine_t G_plain, Pself_plain;
    sm2_jacobian_t self_R;
    uint8_t p_scalar[32], k[32];
    uint8_t self_x[32], self_y[32], gm_x[32], gm_y[32];

    SM2_Z256_POINT gm_P, gm_R;
    sm2_z256_t gp_scalar, gk;

    sm2_get_base_affine(&G_plain);

    memset(p_scalar, 0, sizeof(p_scalar));
    memset(k, 0, sizeof(k));
    fill_pattern_u8(p_scalar, sizeof(p_scalar), seed1);
    fill_pattern_u8(k, sizeof(k), seed2);
    p_scalar[31] |= 1;
    k[31] |= 1;

    /* Local path: P = [p_scalar]G in plain affine, then R = [k]P. */
    sm2_scalar_mul_window_ct_schemeB_mont_to_affine(&Pself_plain, &G_plain, p_scalar);
    sm2_scalar_mul_window_ct_schemeB_mont(&self_R, &Pself_plain, k);

    /* GmSSL path: P = [p_scalar]G, then R = [k]P. */
    sm2_z256_from_bytes(gp_scalar, p_scalar);
    sm2_z256_from_bytes(gk, k);
    sm2_z256_point_mul_generator(&gm_P, gp_scalar);
    sm2_z256_point_mul(&gm_R, gk, &gm_P);

    if (!self_jacobian_mont_to_bytes32(self_x, self_y, &self_R)) {
        printf("self point_mul export failed, seeds=(%u,%u)\n", seed1, seed2);
        return 0;
    }

    if (!gmssl_point_to_bytes32(gm_x, gm_y, &gm_R)) {
        printf("gmssl point_mul export failed, seeds=(%u,%u)\n", seed1, seed2);
        return 0;
    }

    if (!bytes_equal(self_x, gm_x, 32) || !bytes_equal(self_y, gm_y, 32)) {
        printf("Mismatch in point_mul, seeds=(%u,%u)\n", seed1, seed2);
        print_bytes32("self_x", self_x);
        print_bytes32("self_y", self_y);
        print_bytes32("gm_x", gm_x);
        print_bytes32("gm_y", gm_y);
        return 0;
    }

    return 1;
}

static void test_point_mul(void) {
    int ok = 1;
    for (int i = 0; i < TEST_ROUNDS; i++) {
        if (!test_point_mul_once((uint8_t)(i + 17), (uint8_t)(i + 53))) {
            ok = 0;
            break;
        }
    }
    mark_result("point_mul", ok);
}

/* =========================================================
 * 5. scalar stress tests for the scalar-specific mixed-add path
 * ========================================================= */

static int test_point_mul_generator_scalar(const uint8_t d[32], const char *label) {
    sm2_affine_t Gm;
    sm2_jacobian_t self_R;
    SM2_Z256_POINT gm_R;
    sm2_z256_t gd;

    sm2_get_base_affine_mont(&Gm);
    sm2_z256_from_bytes(gd, d);

    sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&self_R, &Gm, d);
    sm2_z256_point_mul_generator(&gm_R, gd);

    return compare_self_gmssl_points(label, &self_R, &gm_R, d, NULL);
}

static int test_point_mul_scalar_pair(
    const uint8_t p_scalar[32],
    const uint8_t k[32],
    const char *label)
{
    sm2_affine_t G_plain, Pself_plain;
    sm2_jacobian_t self_R;
    SM2_Z256_POINT gm_P, gm_R;
    sm2_z256_t gp_scalar, gk;

    sm2_get_base_affine(&G_plain);

    sm2_scalar_mul_window_ct_schemeB_mont_to_affine(&Pself_plain, &G_plain, p_scalar);
    sm2_scalar_mul_window_ct_schemeB_mont(&self_R, &Pself_plain, k);

    sm2_z256_from_bytes(gp_scalar, p_scalar);
    sm2_z256_from_bytes(gk, k);
    sm2_z256_point_mul_generator(&gm_P, gp_scalar);
    sm2_z256_point_mul(&gm_R, gk, &gm_P);

    return compare_self_gmssl_points(label, &self_R, &gm_R, p_scalar, k);
}

static size_t build_stress_scalars(uint8_t scalars[][32], size_t cap) {
    size_t n = 0;

#define ADD_SCALAR(expr) do { \
        if (n < cap) { expr; n++; } \
    } while (0)

    ADD_SCALAR(scalar_set_u64(scalars[n], 0));
    ADD_SCALAR(scalar_set_u64(scalars[n], 1));
    ADD_SCALAR(scalar_set_u64(scalars[n], 2));
    ADD_SCALAR(scalar_set_u64(scalars[n], 3));
    ADD_SCALAR(scalar_set_u64(scalars[n], 4));
    ADD_SCALAR(scalar_set_u64(scalars[n], 5));
    ADD_SCALAR(scalar_set_u64(scalars[n], 16));
    ADD_SCALAR(scalar_set_u64(scalars[n], 17));
    ADD_SCALAR(scalar_set_u64(scalars[n], 31));
    ADD_SCALAR(scalar_set_u64(scalars[n], 32));
    ADD_SCALAR(scalar_set_u64(scalars[n], 33));

    ADD_SCALAR(scalar_sub_small(scalars[n], SM2_N_BYTES, 2));
    ADD_SCALAR(scalar_sub_small(scalars[n], SM2_N_BYTES, 1));
    ADD_SCALAR(scalar_copy(scalars[n], SM2_N_BYTES));
    ADD_SCALAR(scalar_add_small(scalars[n], SM2_N_BYTES, 1));
    ADD_SCALAR(scalar_add_small(scalars[n], SM2_N_BYTES, 2));

    for (int bit = 0; bit < 256; bit += 5) {
        ADD_SCALAR(scalar_set_bit(scalars[n], bit));
    }
    for (int bit = 1; bit < 256; bit += 17) {
        ADD_SCALAR(scalar_set_bit(scalars[n], bit));
    }
    for (int pattern = 0; pattern < 12; pattern++) {
        ADD_SCALAR(scalar_fill_pattern(scalars[n], pattern));
    }

#undef ADD_SCALAR
    return n;
}

static void test_scalar_stress_generator(void) {
    uint8_t scalars[128][32];
    size_t count = build_stress_scalars(scalars, 128);
    int ok = 1;

    for (size_t i = 0; i < count; i++) {
        if (!test_point_mul_generator_scalar(scalars[i], "stress_point_mul_generator")) {
            printf("stress generator index = %zu\n", i);
            ok = 0;
            break;
        }
    }

    for (int i = 0; ok && i < STRESS_ROUNDS; i++) {
        uint8_t d[32];
        fill_pattern_u8(d, 32, (uint8_t)(0x90 + i));
        if ((i & 3) == 0) d[31] &= 0xe0; /* force low 5-bit zero window sometimes */
        if ((i & 7) == 0) d[0] |= 0x80;  /* force high-bit windows */
        if (!test_point_mul_generator_scalar(d, "stress_point_mul_generator_randomish")) {
            printf("stress generator randomish index = %d\n", i);
            ok = 0;
            break;
        }
    }

    mark_result("stress point_mul_generator edge/pattern", ok);
}

static void test_scalar_stress_point_mul(void) {
    uint8_t scalars[128][32];
    size_t count = build_stress_scalars(scalars, 128);
    int ok = 1;

    for (size_t i = 0; ok && i < count; i++) {
        if (scalar_is_zero_bytes(scalars[i]) || scalar_is_order_bytes(scalars[i])) {
            continue; /* GmSSL point_mul input point may not accept infinity as P. */
        }
        for (size_t j = 0; j < count; j += 7) {
            if (!test_point_mul_scalar_pair(scalars[i], scalars[j], "stress_point_mul_pair")) {
                printf("stress point pair indices = (%zu,%zu)\n", i, j);
                ok = 0;
                break;
            }
        }
    }

    for (int i = 0; ok && i < STRESS_ROUNDS; i++) {
        uint8_t p_scalar[32], k[32];
        fill_pattern_u8(p_scalar, 32, (uint8_t)(0x21 + i));
        fill_pattern_u8(k, 32, (uint8_t)(0xc3 - i));
        if ((i & 1) == 0) p_scalar[31] &= 0xe0;
        if ((i & 2) == 0) k[31] &= 0xe0;
        if ((i & 4) == 0) p_scalar[0] |= 0x80;
        if ((i & 8) == 0) k[0] |= 0x80;

        if (!test_point_mul_scalar_pair(p_scalar, k, "stress_point_mul_randomish")) {
            printf("stress point randomish index = %d\n", i);
            ok = 0;
            break;
        }
    }

    mark_result("stress point_mul edge/pattern", ok);
}


/* =========================================================
 * Main entry
 * ========================================================= */

int main(void) {
    printf("Correctness compare with GmSSL\n");
    printf("====================================\n");
    printf("rounds = %d\n\n", TEST_ROUNDS);

    test_fp_to_from_mont();
    test_fp_ops();
    test_fn_ops();
    test_point_mul_generator();
    test_point_mul();
    test_scalar_stress_generator();
    test_scalar_stress_point_mul();

    printf("\n====================================\n");
    printf("TOTAL: %d\n", g_total);
    printf("PASS : %d\n", g_pass);
    printf("FAIL : %d\n", g_fail);

    if (g_fail == 0) {
        printf("[ALL PASS]\n");
        return 0;
    } else {
        printf("[HAS FAILURES]\n");
        return 1;
    }
}
