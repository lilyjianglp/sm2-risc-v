.option push
.option norvc
.text
.align  2
.globl  fp_mont_sqr
.type   fp_mont_sqr, @function


/* ============================================================
 * COND_SUB_P v3（优化版）
 * ============================================================ */
.macro COND_SUB_P r0, r1, r2, r3, c0, c1, c2, c3, borrow, mask

    addi    \c2, \r2, 1

    addi    \c0,     \r0, 1
    sltu    \borrow, \c0, \r0

    sub     \c1,     \r1, a1
    sltu    t6,      \r1, a1
    sub     \c1,     \c1, \borrow
    sltu    \borrow, \c1, \borrow
    or      \borrow, \borrow, t6

    sub     \c2,     \c2, \borrow
    sltu    \borrow, \c2, \borrow

    sub     \c3,     \r3, a2
    sltu    t6,      \r3, a2
    sub     \c3,     \c3, \borrow
    sltu    \borrow, \c3, \borrow
    or      \borrow, \borrow, t6

    xori    \borrow, \borrow, 1
    neg     \mask,   \borrow

    xor     t6, \r0, \c0
    and     t6, t6, \mask
    xor     \r0, \r0, t6

    xor     t6, \r1, \c1
    and     t6, t6, \mask
    xor     \r1, \r1, t6

    xor     t6, \r2, \c2
    and     t6, t6, \mask
    xor     \r2, \r2, t6

    xor     t6, \r3, \c3
    and     t6, t6, \mask
    xor     \r3, \r3, t6

.endm


/* ============================================================
 * ADD 宏
 * ============================================================ */
.macro ADD128_AT_INPLACE l0, l1, l2, top, lo, hi, tmp, c1, c2
    add     \tmp, \l0, \lo
    sltu    \c1,  \tmp, \l0
    mv      \l0,  \tmp

    add     \tmp, \l1, \hi
    sltu    \c2,  \tmp, \l1
    add     \tmp, \tmp, \c1
    sltu    \c1,  \tmp, \c1
    mv      \l1,  \tmp
    or      \c2,  \c2, \c1

    add     \tmp, \l2, \c2
    sltu    \c1,  \tmp, \l2
    mv      \l2,  \tmp

    add     \top, \top, \c1
.endm

.macro ADD128X2_AT_INPLACE l0, l1, l2, top, lo, hi, c1, c2, tmp, topbit
    srli    \topbit, \hi, 63
    srli    \c1,     \lo, 63
    sh1add  \lo,     \lo, zero
    sh1add  \hi,     \hi, \c1

    add     \tmp, \l0, \lo
    sltu    \c1,  \tmp, \l0
    mv      \l0,  \tmp

    add     \tmp, \l1, \hi
    sltu    \c2,  \tmp, \l1
    add     \tmp, \tmp, \c1
    sltu    \c1,  \tmp, \c1
    mv      \l1,  \tmp
    or      \c2,  \c2, \c1

    add     \tmp, \l2, \topbit
    sltu    \c1,  \tmp, \l2
    add     \tmp, \tmp, \c2
    sltu    \c2,  \tmp, \c2
    mv      \l2,  \tmp
    or      \c1,  \c1, \c2

    add     \top, \top, \c1
.endm


/* ============================================================
 * REDC
 * ============================================================ */
.macro REDC_STEP m, x1, x2, x3, x4, x5, top, outreg
    mv      a2, \m
    slli.uw a3, a2, 32
    srli    a4, a2, 32

    mv      \m, zero

    add     a5, \x1, a2
    sltu    a6, a5, \x1
    sltu    a7, a5, a3
    sub     \x1, a5, a3

    add     a5, \x2, a6
    sltu    a6, a5, \x2
    sltu    a0, a5, a4
    sub     a5, a5, a4
    sltu    a1, a5, a7
    sub     \x2, a5, a7
    or      a7, a0, a1

    add     a5, \x3, a6
    sltu    a6, a5, \x3
    sltu    a0, a5, a3
    sub     a5, a5, a3
    sltu    a1, a5, a7
    sub     \x3, a5, a7
    or      a7, a0, a1

    add     a5, \x4, a2
    sltu    a0, a5, \x4
    add     a5, a5, a6
    sltu    a1, a5, a6
    or      a6, a0, a1
    sltu    a0, a5, a4
    sub     a5, a5, a4
    sltu    a1, a5, a7
    sub     \x4, a5, a7
    or      a7, a0, a1

    add     a5, \x5, a6
    sltu    a0, a5, \x5
    sltu    a1, a5, a7
    sub     \x5, a5, a7

    sub     \top, \top, a1
    add     \outreg, \top, a0
.endm


.macro REDC_STEP_LAST m, x1, x2, x3, x4, x5, top, outreg
    REDC_STEP \m, \x1, \x2, \x3, \x4, \x5, \top, \outreg
.endm


/* ============================================================
 * 主函数
 * ============================================================ */
fp_mont_sqr:

    addi    sp, sp, -40
    sd      s0,  0(sp)
    sd      s1,  8(sp)
    sd      s2, 16(sp)
    sd      s3, 24(sp)
    sd      s4, 32(sp)

    mv      s4, a0

    ld      s0,  0(a1)
    ld      s1,  8(a1)
    ld      s2, 16(a1)
    ld      s3, 24(a1)

    mv      t0, zero
    mv      t1, zero
    mv      t2, zero
    mv      t3, zero
    mv      t4, zero
    mv      t5, zero
    mv      t6, zero

    /* ===== round0 ===== */
    mul     a2, s0, s0
    mulhu   a3, s0, s0
    mul     a4, s0, s1
    mulhu   a5, s0, s1
    ADD128_AT_INPLACE t0, t1, t2, t6, a2, a3, a6, a7, a0

    mul     a0, s0, s2
    mulhu   a1, s0, s2
    ADD128X2_AT_INPLACE t1, t2, t3, t6, a4, a5, a6, a7, a2, a3

    mul     a4, s0, s3
    mulhu   a5, s0, s3
    ADD128X2_AT_INPLACE t2, t3, t4, t6, a0, a1, a6, a7, a2, a3
    ADD128X2_AT_INPLACE t3, t4, t5, t6, a4, a5, a6, a7, a2, a3

    REDC_STEP t0, t1, t2, t3, t4, t5, t6, t0
    mv t6, zero

    /* ===== round1 ===== */
    mul     a2, s1, s1
    mulhu   a3, s1, s1
    mul     a4, s1, s2
    mulhu   a5, s1, s2
    ADD128_AT_INPLACE t2, t3, t4, t6, a2, a3, a6, a7, a0

    mul     a0, s1, s3
    mulhu   a1, s1, s3
    ADD128X2_AT_INPLACE t3, t4, t5, t6, a4, a5, a6, a7, a2, a3
    ADD128X2_AT_INPLACE t4, t5, t0, t6, a0, a1, a6, a7, a2, a3

    REDC_STEP t1, t2, t3, t4, t5, t0, t6, t1
    mv t6, zero

    /* ===== round2 ===== */
    mul     a2, s2, s2
    mulhu   a3, s2, s2
    mul     a4, s2, s3
    mulhu   a5, s2, s3
    ADD128_AT_INPLACE t4, t5, t0, t6, a2, a3, a6, a7, a0
    ADD128X2_AT_INPLACE t5, t0, t1, t6, a4, a5, a6, a7, a0, a1

    REDC_STEP t2, t3, t4, t5, t0, t1, t6, t2
    mv t6, zero

    /* ===== round3 ===== */
    mul     a2, s3, s3
    mulhu   a3, s3, s3
    ADD128_AT_INPLACE t0, t1, t2, t6, a2, a3, a6, a7, a0

    REDC_STEP_LAST t3, t4, t5, t0, t1, t2, t6, t3

    /* ===== tail ===== */
    add     a7, t2, t3
    li      a0, -1
    andi    a6, a7, 1
    slli    a1, a0, 32
    slli    a3, a6, 32
    addi    a2, a1, -1
    neg     t2, a6

    add     t4, t4, a6
    sltu    t3, t4, a6

    srli    a5, t2, 32
    add     a4, t5, a5
    sltu    a7, a4, t5
    add     t5, a4, t3
    sltu    a5, t5, t3
    or      t3, a7, a5

    add     t0, t0, t3
    sltu    t3, t0, t3

    add     a4, t1, a3
    sltu    a7, a4, t1
    add     t1, a4, t3
    sltu    a5, t1, t3
    or      t3, a7, a5

    COND_SUB_P t4, t5, t0, t1, a3, a4, a5, a6, t2, t3

    sd      t4,  0(s4)
    sd      t5,  8(s4)
    sd      t0, 16(s4)
    sd      t1, 24(s4)

    ld      s0,  0(sp)
    ld      s1,  8(sp)
    ld      s2, 16(sp)
    ld      s3, 24(sp)
    ld      s4, 32(sp)
    addi    sp, sp, 40
    ret

.size fp_mont_sqr, .-fp_mont_sqr
.option pop