    .text
    .align 2

/*
 * Conditional subtraction by SM2 p.
 *
 * Input/output:
 *   s0,s1,s2,s3 = 256-bit value, little-endian limbs
 *
 * Modulus:
 *   t0 = p0 = 0xFFFFFFFFFFFFFFFF
 *   t1 = p1 = 0xFFFFFFFF00000000
 *   t2 = p2 = 0xFFFFFFFFFFFFFFFF
 *   t3 = p3 = 0xFFFFFFFEFFFFFFFF
 *
 * Clobbers:
 *   a0,a1,a2,a3,a4,a5,a6,t4
 *
 * Correct borrow logic:
 *   d  = a - b
 *   b1 = (a < b)
 *   b2 = (d < borrow)     // must be computed before d -= borrow
 *   d -= borrow
 *   borrow = b1 | b2
 */
    .macro COND_SUB_SM2_P
    /* limb 0: no incoming borrow */
    sub     a0, s0, t0
    sltu    a4, s0, t0        /* borrow */

    /* limb 1 */
    sub     a1, s1, t1        /* d = s1 - p1 */
    sltu    a5, s1, t1        /* b1 */
    sltu    a6, a1, a4        /* b2 = d < borrow, before subtracting borrow */
    sub     a1, a1, a4        /* d -= borrow */
    or      a4, a5, a6        /* new borrow */

    /* limb 2 */
    sub     a2, s2, t2
    sltu    a5, s2, t2
    sltu    a6, a2, a4
    sub     a2, a2, a4
    or      a4, a5, a6

    /* limb 3 */
    sub     a3, s3, t3
    sltu    a5, s3, t3
    sltu    a6, a3, a4
    sub     a3, a3, a4
    or      a4, a5, a6        /* final borrow */

    /*
     * If final borrow == 0, keep d = s - p.
     * If final borrow == 1, keep original s.
     */
    xori    a4, a4, 1
    sub     a4, x0, a4        /* mask = no_borrow ? -1 : 0 */
    xori    t4, a4, -1        /* inverse mask */

    and     a0, a0, a4
    and     a1, a1, a4
    and     a2, a2, a4
    and     a3, a3, a4

    and     s0, s0, t4
    and     s1, s1, t4
    and     s2, s2, t4
    and     s3, s3, t4

    or      s0, s0, a0
    or      s1, s1, a1
    or      s2, s2, a2
    or      s3, s3, a3
    .endm

/*
 * Store r = -a mod p into dst.
 *
 * Input:
 *   a0 = dst
 *   a1 = a
 *
 * Clobbers:
 *   a2,a3,a4,a5,a6,a7,t0,t1,t2,t3,t4,t5,t6
 */
    .macro NEG_MOD_SM2 dst, src
    li      t0, -1
    li      t1, 0xFFFFFFFF00000000
    li      t2, -1
    li      t3, 0xFFFFFFFEFFFFFFFF

    ld      a2, 0(\src)
    ld      a3, 8(\src)
    ld      a4, 16(\src)
    ld      a5, 24(\src)

    or      t4, a2, a3
    or      t4, t4, a4
    or      t4, t4, a5
    sltu    t4, x0, t4
    sub     t4, x0, t4       /* nonzero mask */

    sub     t5, t0, a2
    sltu    t6, t0, a2       /* borrow */

    sub     a6, t1, a3
    sltu    a7, t1, a3
    sltu    a0, a6, t6
    sub     a6, a6, t6
    or      t6, a7, a0

    sub     a7, t2, a4
    sltu    a0, t2, a4
    sltu    a1, a7, t6
    sub     a7, a7, t6
    or      t6, a0, a1

    sub     a0, t3, a5
    sltu    a1, t3, a5
    sltu    a2, a0, t6
    sub     a0, a0, t6

    and     t5, t5, t4
    and     a6, a6, t4
    and     a7, a7, t4
    and     a0, a0, t4

    sd      t5, 0(\dst)
    sd      a6, 8(\dst)
    sd      a7, 16(\dst)
    sd      a0, 24(\dst)
    .endm

/*
 * Expand one:
 *   r = c1*x + c2*y mod p
 *
 * Arguments are stack offsets in the current frame:
 *   dstoff, c1off, xoff, negxoff, c2off, yoff, negyoff
 *
 * Scratch:
 *   uses t[5] at 64(sp), clobbers s0-s6, a0-a7, t0-t6.
 */
    .macro LINCOMB_FROM_STACK dstoff, c1off, xoff, negxoff, c2off, yoff, negyoff
    ld      s6, \dstoff(sp)
    ld      a1, \c1off(sp)
    ld      a2, \xoff(sp)
    addi    a3, sp, \negxoff
    ld      a4, \c2off(sp)
    ld      a5, \yoff(sp)
    addi    a6, sp, \negyoff
    addi    s5, sp, 64

    /* Select source pointers and absolute coefficients. */
    srai    t0, a1, 63
    xor     t1, a2, a3
    and     t1, t1, t0
    xor     s0, a2, t1
    xor     s2, a1, t0
    sub     s2, s2, t0

    srai    t0, a4, 63
    xor     t1, a5, a6
    and     t1, t1, t0
    xor     s1, a5, t1
    xor     s3, a4, t0
    sub     s3, s3, t0

    li      s4, 0

    /* limb 0 */
    ld      t0, 0(s0)
    ld      t1, 0(s1)
    mul     t2, s2, t0
    mulhu   t3, s2, t0
    mul     t4, s3, t1
    mulhu   t5, s3, t1
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     t6, t6, s4
    sltu    a1, t6, s4
    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1
    sd      t6, 0(s5)

    /* limb 1 */
    ld      t0, 8(s0)
    ld      t1, 8(s1)
    mul     t2, s2, t0
    mulhu   t3, s2, t0
    mul     t4, s3, t1
    mulhu   t5, s3, t1
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     t6, t6, s4
    sltu    a1, t6, s4
    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1
    sd      t6, 8(s5)

    /* limb 2 */
    ld      t0, 16(s0)
    ld      t1, 16(s1)
    mul     t2, s2, t0
    mulhu   t3, s2, t0
    mul     t4, s3, t1
    mulhu   t5, s3, t1
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     t6, t6, s4
    sltu    a1, t6, s4
    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1
    sd      t6, 16(s5)

    /* limb 3 */
    ld      t0, 24(s0)
    ld      t1, 24(s1)
    mul     t2, s2, t0
    mulhu   t3, s2, t0
    mul     t4, s3, t1
    mulhu   t5, s3, t1
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     t6, t6, s4
    sltu    a1, t6, s4
    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1
    sd      t6, 24(s5)
    sd      s4, 32(s5)

    /* Load t into fold state. */
    ld      s0, 0(s5)
    ld      s1, 8(s5)
    ld      s2, 16(s5)
    ld      s3, 24(s5)
    ld      s4, 32(s5)

    /* fold round 1 */
    add     t0, s0, s4
    sltu    t1, t0, s0
    mv      s0, t0
    slli    t2, s4, 32
    srli    t3, s4, 32
    sub     t4, t2, s4
    sltu    t5, t2, s4
    sub     t6, t3, t5
    add     a0, s1, t4
    sltu    a1, a0, s1
    add     a0, a0, t1
    sltu    a2, a0, t1
    mv      s1, a0
    add     a3, t6, a1
    add     a3, a3, a2
    add     a4, s2, t3
    sltu    a5, a4, s2
    add     a4, a4, a3
    sltu    a6, a4, a3
    mv      s2, a4
    add     a7, a5, a6
    slli    t2, s4, 32
    srli    t3, s4, 32
    add     a0, s3, t2
    sltu    a1, a0, s3
    add     a0, a0, a7
    sltu    a2, a0, a7
    mv      s3, a0
    add     s4, t3, a1
    add     s4, s4, a2

    /* fold round 2 */
    add     t0, s0, s4
    sltu    t1, t0, s0
    mv      s0, t0
    slli    t2, s4, 32
    srli    t3, s4, 32
    sub     t4, t2, s4
    sltu    t5, t2, s4
    sub     t6, t3, t5
    add     a0, s1, t4
    sltu    a1, a0, s1
    add     a0, a0, t1
    sltu    a2, a0, t1
    mv      s1, a0
    add     a3, t6, a1
    add     a3, a3, a2
    add     a4, s2, t3
    sltu    a5, a4, s2
    add     a4, a4, a3
    sltu    a6, a4, a3
    mv      s2, a4
    add     a7, a5, a6
    slli    t2, s4, 32
    srli    t3, s4, 32
    add     a0, s3, t2
    sltu    a1, a0, s3
    add     a0, a0, a7
    sltu    a2, a0, a7
    mv      s3, a0
    add     s4, t3, a1
    add     s4, s4, a2

    li      t0, -1
    li      t1, 0xFFFFFFFF00000000
    li      t2, -1
    li      t3, 0xFFFFFFFEFFFFFFFF
    COND_SUB_SM2_P
    COND_SUB_SM2_P
    COND_SUB_SM2_P

    sd      s0, 0(s6)
    sd      s1, 8(s6)
    sd      s2, 16(s6)
    sd      s3, 24(s6)
    .endm

#ifndef PORNIN_MACRO_ONLY
    .globl pornin_update_uv_asm
    .type  pornin_update_uv_asm, @function

/*
 * void pornin_update_uv_asm(
 *     fp_t *new_u, fp_t *new_v,
 *     const fp_t *u, const fp_t *v,
 *     int64_t f0, int64_t g0,
 *     int64_t f1, int64_t g1);
 */
pornin_update_uv_asm:
    /*
     * Fused leaf version. Stack layout:
     *   0..31    neg_u
     *   32..63   neg_v
     *   64..103  t[5] scratch for LINCOMB_FROM_STACK
     *   104..152 saved s0..s6
     *   160      new_u
     *   168      new_v
     *   176      u
     *   184      v
     *   192      f0
     *   200      g0
     *   208      f1
     *   216      g1
     */
    addi    sp, sp, -224

    sd      s0, 104(sp)
    sd      s1, 112(sp)
    sd      s2, 120(sp)
    sd      s3, 128(sp)
    sd      s4, 136(sp)
    sd      s5, 144(sp)
    sd      s6, 152(sp)

    sd      a0, 160(sp)       /* new_u */
    sd      a1, 168(sp)       /* new_v */
    sd      a2, 176(sp)       /* u */
    sd      a3, 184(sp)       /* v */
    sd      a4, 192(sp)       /* f0 */
    sd      a5, 200(sp)       /* g0 */
    sd      a6, 208(sp)       /* f1 */
    sd      a7, 216(sp)       /* g1 */

    addi    s0, sp, 0         /* neg_u */
    ld      s1, 176(sp)       /* u */
    NEG_MOD_SM2 s0, s1

    addi    s0, sp, 32        /* neg_v */
    ld      s1, 184(sp)       /* v */
    NEG_MOD_SM2 s0, s1

    /* new_u = f0*u + g0*v mod p */
    LINCOMB_FROM_STACK 160, 192, 176, 0, 200, 184, 32

    /* new_v = f1*u + g1*v mod p */
    LINCOMB_FROM_STACK 168, 208, 176, 0, 216, 184, 32

    ld      s0, 104(sp)
    ld      s1, 112(sp)
    ld      s2, 120(sp)
    ld      s3, 128(sp)
    ld      s4, 136(sp)
    ld      s5, 144(sp)
    ld      s6, 152(sp)
    addi    sp, sp, 224
    ret

    .size pornin_update_uv_asm, .-pornin_update_uv_asm

    .globl pornin_fp_lincomb_small_mod_precomp_asm
    .type  pornin_fp_lincomb_small_mod_precomp_asm, @function

/*
 * void pornin_fp_lincomb_small_mod_precomp_asm(
 *     fp_t *r,
 *     int64_t c1, const fp_t *x, const fp_t *neg_x,
 *     int64_t c2, const fp_t *y, const fp_t *neg_y);
 *
 * ABI:
 *   a0 = r
 *   a1 = c1
 *   a2 = x
 *   a3 = neg_x
 *   a4 = c2
 *   a5 = y
 *   a6 = neg_y
 *
 * Semantics:
 *   r = c1*x + c2*y mod p
 *
 * Method:
 *   1. Select x or neg_x according to sign(c1)
 *   2. Select y or neg_y according to sign(c2)
 *   3. Compute t[0..4] = abs(c1)*sx + abs(c2)*sy
 *   4. Fold 320-bit t modulo SM2 p
 *   5. Apply three fixed conditional subtractions by p
 */
pornin_fp_lincomb_small_mod_precomp_asm:
    addi    sp, sp, -112

    sd      ra, 40(sp)
    sd      s0, 48(sp)
    sd      s1, 56(sp)
    sd      s2, 64(sp)
    sd      s3, 72(sp)
    sd      s4, 80(sp)
    sd      s5, 88(sp)
    sd      s6, 96(sp)

    /* s6 = output pointer r, s5 = local t[5] base */
    mv      s6, a0
    mv      s5, sp

    /*
     * Select source pointers in constant time.
     *
     * s0 = (c1 < 0) ? neg_x : x
     * s1 = (c2 < 0) ? neg_y : y
     * s2 = abs(c1)
     * s3 = abs(c2)
     */
    srai    t0, a1, 63         /* mask1 = c1 < 0 ? -1 : 0 */
    xor     t1, a2, a3
    and     t1, t1, t0
    xor     s0, a2, t1         /* selected x pointer */

    xor     s2, a1, t0
    sub     s2, s2, t0         /* abs(c1) */

    srai    t0, a4, 63         /* mask2 = c2 < 0 ? -1 : 0 */
    xor     t1, a5, a6
    and     t1, t1, t0
    xor     s1, a5, t1         /* selected y pointer */

    xor     s3, a4, t0
    sub     s3, s3, t0         /* abs(c2) */

    li      s4, 0              /* carry */

    /*
     * Compute:
     *   t = abs(c1)*selected_x + abs(c2)*selected_y
     *
     * Since coefficients are already absolute values, use unsigned high multiply.
     */

    /* limb 0 */
    ld      t0, 0(s0)
    ld      t1, 0(s1)

    mul     t2, s2, t0
    mulhu   t3, s2, t0

    mul     t4, s3, t1
    mulhu   t5, s3, t1

    add     t6, t2, t4
    sltu    a0, t6, t2

    add     t6, t6, s4
    sltu    a1, t6, s4

    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1

    sd      t6, 0(s5)

    /* limb 1 */
    ld      t0, 8(s0)
    ld      t1, 8(s1)

    mul     t2, s2, t0
    mulhu   t3, s2, t0

    mul     t4, s3, t1
    mulhu   t5, s3, t1

    add     t6, t2, t4
    sltu    a0, t6, t2

    add     t6, t6, s4
    sltu    a1, t6, s4

    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1

    sd      t6, 8(s5)

    /* limb 2 */
    ld      t0, 16(s0)
    ld      t1, 16(s1)

    mul     t2, s2, t0
    mulhu   t3, s2, t0

    mul     t4, s3, t1
    mulhu   t5, s3, t1

    add     t6, t2, t4
    sltu    a0, t6, t2

    add     t6, t6, s4
    sltu    a1, t6, s4

    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1

    sd      t6, 16(s5)

    /* limb 3 */
    ld      t0, 24(s0)
    ld      t1, 24(s1)

    mul     t2, s2, t0
    mulhu   t3, s2, t0

    mul     t4, s3, t1
    mulhu   t5, s3, t1

    add     t6, t2, t4
    sltu    a0, t6, t2

    add     t6, t6, s4
    sltu    a1, t6, s4

    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1

    sd      t6, 24(s5)
    sd      s4, 32(s5)

    /*
     * Inline pornin_sm2_fold_u320().
     *
     * State:
     *   s0,s1,s2,s3 = x0,x1,x2,x3
     *   s4          = h
     *
     * SM2 relation:
     *   2^256 = 2^224 + 2^96 - 2^64 + 1 mod p
     */
    ld      s0, 0(s5)
    ld      s1, 8(s5)
    ld      s2, 16(s5)
    ld      s3, 24(s5)
    ld      s4, 32(s5)

    /* ---------------- fold round 1 ---------------- */

    /* x0 += h */
    add     t0, s0, s4
    sltu    t1, t0, s0        /* carry from x0 */
    mv      s0, t0

    /*
     * expr = ((uint128)h << 32) - h
     * low(expr)  = t4
     * high(expr) = t6
     */
    slli    t2, s4, 32
    srli    t3, s4, 32
    sub     t4, t2, s4
    sltu    t5, t2, s4
    sub     t6, t3, t5

    /* x1 += low(expr) + carry_x0 */
    add     a0, s1, t4
    sltu    a1, a0, s1

    add     a0, a0, t1
    sltu    a2, a0, t1

    mv      s1, a0

    /* carry into x2 from x1 plus high(expr) */
    add     a3, t6, a1
    add     a3, a3, a2

    /* x2 += (h >> 32) + carry */
    add     a4, s2, t3
    sltu    a5, a4, s2

    add     a4, a4, a3
    sltu    a6, a4, a3

    mv      s2, a4
    add     a7, a5, a6        /* carry into x3 */

    /* x3 += (h << 32) + carry */
    slli    t2, s4, 32
    srli    t3, s4, 32

    add     a0, s3, t2
    sltu    a1, a0, s3

    add     a0, a0, a7
    sltu    a2, a0, a7

    mv      s3, a0

    /* new h = (h >> 32) + carry_x3 */
    add     s4, t3, a1
    add     s4, s4, a2

    /* ---------------- fold round 2 ---------------- */

    /* x0 += h */
    add     t0, s0, s4
    sltu    t1, t0, s0
    mv      s0, t0

    /*
     * expr = ((uint128)h << 32) - h
     */
    slli    t2, s4, 32
    srli    t3, s4, 32
    sub     t4, t2, s4
    sltu    t5, t2, s4
    sub     t6, t3, t5

    /* x1 += low(expr) + carry_x0 */
    add     a0, s1, t4
    sltu    a1, a0, s1

    add     a0, a0, t1
    sltu    a2, a0, t1

    mv      s1, a0

    /* carry into x2 from x1 plus high(expr) */
    add     a3, t6, a1
    add     a3, a3, a2

    /* x2 += (h >> 32) + carry */
    add     a4, s2, t3
    sltu    a5, a4, s2

    add     a4, a4, a3
    sltu    a6, a4, a3

    mv      s2, a4
    add     a7, a5, a6

    /* x3 += (h << 32) + carry */
    slli    t2, s4, 32
    srli    t3, s4, 32

    add     a0, s3, t2
    sltu    a1, a0, s3

    add     a0, a0, a7
    sltu    a2, a0, a7

    mv      s3, a0

    /* new h = (h >> 32) + carry_x3 */
    add     s4, t3, a1
    add     s4, s4, a2

    /*
     * Three fixed conditional subtractions by p.
     *
     * SM2 p limbs, little-endian:
     *   p0 = 0xFFFFFFFFFFFFFFFF
     *   p1 = 0xFFFFFFFF00000000
     *   p2 = 0xFFFFFFFFFFFFFFFF
     *   p3 = 0xFFFFFFFEFFFFFFFF
     */
    li      t0, -1
    li      t1, 0xFFFFFFFF00000000
    li      t2, -1
    li      t3, 0xFFFFFFFEFFFFFFFF

    COND_SUB_SM2_P
    COND_SUB_SM2_P
    COND_SUB_SM2_P

    /* Store result */
    sd      s0, 0(s6)
    sd      s1, 8(s6)
    sd      s2, 16(s6)
    sd      s3, 24(s6)

    /* Restore callee-saved registers */
    ld      ra, 40(sp)
    ld      s0, 48(sp)
    ld      s1, 56(sp)
    ld      s2, 64(sp)
    ld      s3, 72(sp)
    ld      s4, 80(sp)
    ld      s5, 88(sp)
    ld      s6, 96(sp)

    addi    sp, sp, 112
    ret

    .size pornin_fp_lincomb_small_mod_precomp_asm, .-pornin_fp_lincomb_small_mod_precomp_asm
#endif
