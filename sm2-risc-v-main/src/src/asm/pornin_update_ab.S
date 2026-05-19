    .text
    .align 2

/*
 * pornin_update_ab.S v3
 *
 * This version intentionally uses s-registers inside the fused update_ab
 * path.  The previous caller-saved-only version had too much register
 * pressure, which made the carry chain very hard to schedule on X60.
 *
 * Exports:
 *   uint64_t pornin_lincomb_shift31_abs_asm(...)
 *   void     pornin_update_ab_asm(...)
 */

/*
 * Input:
 *   a0 = out
 *   a1 = c1          signed int64
 *   a2 = x           const fp_t*
 *   a3 = c2          signed int64
 *   a4 = y           const fp_t*
 *
 * Output:
 *   a0 = is_neg
 *
 * Clobbers:
 *   a1-a7, t0-t6, s0-s5
 *
 * Register plan:
 *   s0 = out
 *   s1 = z0
 *   s2 = z1
 *   s3 = z2
 *   s4 = z3
 *   s5 = z4 / signed carry
 */
    .macro LINCOMB_SHIFT31_ABS_SREG
    mv      s0, a0

    /* ---------- limb 0 ---------- */
    ld      t0, 0(a2)
    ld      t1, 0(a4)
    mul     t2, a1, t0
    mulhsu  t3, a1, t0
    mul     t4, a3, t1
    mulhsu  t5, a3, t1

    ld      a6, 8(a2)
    ld      a7, 8(a4)

    add     t6, t2, t4
    sltu    a5, t6, t2
    add     s5, t3, t5
    add     s5, s5, a5
    mv      s1, t6

    /* ---------- limb 1 ---------- */
    mul     t2, a1, a6
    mulhsu  t3, a1, a6
    mul     t4, a3, a7
    mulhsu  t5, a3, a7

    add     t6, t2, t4
    sltu    a5, t6, t2
    add     t6, t6, s5
    sltu    a6, t6, s5
    srai    a7, s5, 63

    add     s5, t3, t5
    add     s5, s5, a5
    add     s5, s5, a6
    add     s5, s5, a7
    mv      s2, t6

    /* ---------- limb 2 ---------- */
    ld      t0, 16(a2)
    ld      t1, 16(a4)
    mul     t2, a1, t0
    mulhsu  t3, a1, t0
    mul     t4, a3, t1
    mulhsu  t5, a3, t1

    add     t6, t2, t4
    sltu    a5, t6, t2
    add     t6, t6, s5
    sltu    a6, t6, s5
    srai    a7, s5, 63

    add     s5, t3, t5
    add     s5, s5, a5
    add     s5, s5, a6
    add     s5, s5, a7
    mv      s3, t6

    /* ---------- limb 3 ---------- */
    ld      t0, 24(a2)
    ld      t1, 24(a4)
    mul     t2, a1, t0
    mulhsu  t3, a1, t0
    mul     t4, a3, t1
    mulhsu  t5, a3, t1

    add     t6, t2, t4
    sltu    a5, t6, t2
    add     t6, t6, s5
    sltu    a6, t6, s5
    srai    a7, s5, 63

    add     s5, t3, t5
    add     s5, s5, a5
    add     s5, s5, a6
    add     s5, s5, a7
    mv      s4, t6

    /*
     * Absolute value over five limbs s1,s2,s3,s4,s5.
     */
    srli    a7, s5, 63         /* is_neg */
    sub     a6, x0, a7         /* mask = -is_neg */
    mv      a1, a7             /* add_one = is_neg */

    xor     s1, s1, a6
    add     s1, s1, a1
    sltu    a1, s1, a1

    xor     s2, s2, a6
    add     s2, s2, a1
    sltu    a1, s2, a1

    xor     s3, s3, a6
    add     s3, s3, a1
    sltu    a1, s3, a1

    xor     s4, s4, a6
    add     s4, s4, a1
    sltu    a1, s4, a1

    xor     s5, s5, a6
    add     s5, s5, a1

    /*
     * Exact right shift by 31.
     */
    srli    a1, s1, 31
    slli    a2, s2, 33
    or      a1, a1, a2
    sd      a1, 0(s0)

    srli    a1, s2, 31
    slli    a2, s3, 33
    or      a1, a1, a2
    sd      a1, 8(s0)

    srli    a1, s3, 31
    slli    a2, s4, 33
    or      a1, a1, a2
    sd      a1, 16(s0)

    srli    a1, s4, 31
    slli    a2, s5, 33
    or      a1, a1, a2
    sd      a1, 24(s0)

    mv      a0, a7
    .endm

    .macro COND_NEG_COEFF_REG ptr
    ld      t2, 0(\ptr)
    sub     t3, x0, t2
    xor     t4, t2, t3
    and     t4, t4, a7
    xor     t2, t2, t4
    sd      t2, 0(\ptr)
    .endm

    .globl pornin_lincomb_shift31_abs_asm
    .type  pornin_lincomb_shift31_abs_asm, @function

pornin_lincomb_shift31_abs_asm:
    addi    sp, sp, -48
    sd      s0, 0(sp)
    sd      s1, 8(sp)
    sd      s2, 16(sp)
    sd      s3, 24(sp)
    sd      s4, 32(sp)
    sd      s5, 40(sp)

    LINCOMB_SHIFT31_ABS_SREG

    ld      s0, 0(sp)
    ld      s1, 8(sp)
    ld      s2, 16(sp)
    ld      s3, 24(sp)
    ld      s4, 32(sp)
    ld      s5, 40(sp)
    addi    sp, sp, 48
    ret

    .size pornin_lincomb_shift31_abs_asm, .-pornin_lincomb_shift31_abs_asm

    .globl pornin_update_ab_asm
    .type  pornin_update_ab_asm, @function

/*
 * void pornin_update_ab_asm(
 *     fp_t *new_a, fp_t *new_b,
 *     const fp_t *a, const fp_t *b,
 *     int64_t *f0, int64_t *g0,
 *     int64_t *f1, int64_t *g1);
 */
pornin_update_ab_asm:
    /*
     * Stack layout, 128 bytes:
     *   0..95    saved s0..s11
     *   96       f1 ptr
     *   104      g1 ptr
     */
    addi    sp, sp, -128
    sd      s0, 0(sp)
    sd      s1, 8(sp)
    sd      s2, 16(sp)
    sd      s3, 24(sp)
    sd      s4, 32(sp)
    sd      s5, 40(sp)
    sd      s6, 48(sp)
    sd      s7, 56(sp)
    sd      s8, 64(sp)
    sd      s9, 72(sp)
    sd      s10, 80(sp)
    sd      s11, 88(sp)
    sd      a6, 96(sp)
    sd      a7, 104(sp)

    mv      s6, a0             /* new_a */
    mv      s7, a1             /* new_b */
    mv      s8, a2             /* a */
    mv      s9, a3             /* b */
    mv      s10, a4            /* f0 ptr */
    mv      s11, a5            /* g0 ptr */

    /*
     * new_a = abs((f0*a + g0*b) >> 31)
     */
    mv      a0, s6
    ld      a1, 0(s10)
    mv      a2, s8
    ld      a3, 0(s11)
    mv      a4, s9
    LINCOMB_SHIFT31_ABS_SREG

    sub     a7, x0, a0
    COND_NEG_COEFF_REG s10
    COND_NEG_COEFF_REG s11

    /*
     * new_b = abs((f1*a + g1*b) >> 31)
     */
    ld      t0, 96(sp)         /* f1 ptr */
    ld      t1, 104(sp)        /* g1 ptr */
    mv      a0, s7
    ld      a1, 0(t0)
    mv      a2, s8
    ld      a3, 0(t1)
    mv      a4, s9
    LINCOMB_SHIFT31_ABS_SREG

    sub     a7, x0, a0
    ld      t0, 96(sp)
    ld      t1, 104(sp)
    COND_NEG_COEFF_REG t0
    COND_NEG_COEFF_REG t1

    ld      s0, 0(sp)
    ld      s1, 8(sp)
    ld      s2, 16(sp)
    ld      s3, 24(sp)
    ld      s4, 32(sp)
    ld      s5, 40(sp)
    ld      s6, 48(sp)
    ld      s7, 56(sp)
    ld      s8, 64(sp)
    ld      s9, 72(sp)
    ld      s10, 80(sp)
    ld      s11, 88(sp)
    addi    sp, sp, 128
    ret

    .size pornin_update_ab_asm, .-pornin_update_ab_asm
