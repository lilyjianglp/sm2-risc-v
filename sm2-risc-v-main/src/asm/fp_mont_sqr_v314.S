    .option push
    .option norvc
    .text
    .align  2
    .globl  fp_mont_sqr
    .type   fp_mont_sqr, @function

/* ============================================================
 * fp_mont_sqr(fp_t* r, const fp_t* a)
 *
 * Diagnostic version:
 * - keep REDC_STEP refactor structure
 * - keep compensation path:
 *     extra+t6 -> mask -> add mask * FP_MONT_ONE
 * - comment out ONLY the final COND_SUB_P
 *
 * WARNING:
 *   This version is for cycle diagnosis only.
 *   Do NOT use it for correctness tests.
 * ============================================================ */

    /* a0 = -1
     * a1 = 0xFFFFFFFF00000000
     * a2 = 0xFFFFFFFEFFFFFFFF
     */
    .macro COND_SUB_P r0, r1, r2, r3, c0, c1, c2, c3, borrow, mask
    mv      \borrow, zero

    sub     t6, \r0, a0
    sltu    a7, \r0, a0
    sub     \c0, t6, \borrow
    sltu    \borrow, t6, \borrow
    or      \borrow, \borrow, a7

    sub     t6, \r1, a1
    sltu    a7, \r1, a1
    sub     \c1, t6, \borrow
    sltu    \borrow, t6, \borrow
    or      \borrow, \borrow, a7

    sub     t6, \r2, a0
    sltu    a7, \r2, a0
    sub     \c2, t6, \borrow
    sltu    \borrow, t6, \borrow
    or      \borrow, \borrow, a7

    sub     t6, \r3, a2
    sltu    a7, \r3, a2
    sub     \c3, t6, \borrow
    sltu    \borrow, t6, \borrow
    or      \borrow, \borrow, a7

    xori    \borrow, \borrow, 1
    neg     \mask, \borrow

    and     \c0, \c0, \mask
    andn    t6,  \r0, \mask
    or      \r0, t6, \c0

    and     \c1, \c1, \mask
    andn    t6,  \r1, \mask
    or      \r1, t6, \c1

    and     \c2, \c2, \mask
    andn    t6,  \r2, \mask
    or      \r2, t6, \c2

    and     \c3, \c3, \mask
    andn    t6,  \r3, \mask
    or      \r3, t6, \c3
    .endm

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
    slli    \lo,     \lo, 1
    slli    \tmp,    \hi, 1
    or      \hi,     \tmp, \c1

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

    /* unified Montgomery reduction step:
     * m  = current low limb to eliminate
     * x1..x5 = next 5 active limbs in circular order
     * top = top carry/borrow accumulator
     * outreg = emitted reduced limb
     *
     * For rounds 0..2, caller clears top afterward.
     * For round 3, caller keeps top because tail still uses it.
     */
       .macro REDC_STEP_CLR m, x1, x2, x3, x4, x5, top, outreg
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
    mv      \top, zero
    .endm

    .macro REDC_STEP_KEEP m, x1, x2, x3, x4, x5, top, outreg
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
fp_mont_sqr:
    addi    sp, sp, -40
    sd      s0,  0(sp)
    sd      s1,  8(sp)
    sd      s2, 16(sp)
    sd      s3, 24(sp)
    sd      a0, 32(sp)          /* spill output pointer */

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

    /* round 0 */
    mul     a2, s0, s0
    mulhu   a3, s0, s0
    mul     a4, s0, s1
    mulhu   a5, s0, s1
    ADD128_AT_INPLACE   t0, t1, t2, t6, a2, a3, a6, a7, a0

    mul     a0, s0, s2
    mulhu   a1, s0, s2
    ADD128X2_AT_INPLACE t1, t2, t3, t6, a4, a5, a6, a7, a2, a3

    mul     a4, s0, s3
    mulhu   a5, s0, s3
    ADD128X2_AT_INPLACE t2, t3, t4, t6, a0, a1, a6, a7, a2, a3
    ADD128X2_AT_INPLACE t3, t4, t5, t6, a4, a5, a6, a7, a2, a3

    REDC_STEP_CLR t0, t1, t2, t3, t4, t5, t6, t0

    /* round 1 */
    mul     a2, s1, s1
    mulhu   a3, s1, s1
    mul     a4, s1, s2
    mulhu   a5, s1, s2
    ADD128_AT_INPLACE   t2, t3, t4, t6, a2, a3, a6, a7, a0

    mul     a0, s1, s3
    mulhu   a1, s1, s3
    ADD128X2_AT_INPLACE t3, t4, t5, t6, a4, a5, a6, a7, a2, a3
    ADD128X2_AT_INPLACE t4, t5, t0, t6, a0, a1, a6, a7, a2, a3

    REDC_STEP_CLR t1, t2, t3, t4, t5, t0, t6, t1

    /* round 2 */
    mul     a2, s2, s2
    mulhu   a3, s2, s2
    mul     a4, s2, s3
    mulhu   a5, s2, s3
    ADD128_AT_INPLACE   t4, t5, t0, t6, a2, a3, a6, a7, a0
    ADD128X2_AT_INPLACE t5, t0, t1, t6, a4, a5, a6, a7, a0, a1

    REDC_STEP_CLR t2, t3, t4, t5, t0, t1, t6, t2

    /* round 3 */
    mul     a2, s3, s3
    mulhu   a3, s3, s3
    ADD128_AT_INPLACE   t0, t1, t2, t6, a2, a3, a6, a7, a0

    REDC_STEP_KEEP t3, t4, t5, t0, t1, t2, t6, t6

    /* final logical result:
     * discard w0=t3
     * r0=t4, r1=t5, r2=t0, r3=t1
     * extra=t2+t6
     */
    add     t2, t2, t6

    /* cond = extra & 1, mask = -cond */
    andi    t2, t2, 1
    neg     t2, t2

    /* preload FP_P limbs early for final subtract
     * a0 = 0xFFFFFFFFFFFFFFFF
     * a1 = 0xFFFFFFFF00000000
     * a2 = 0xFFFFFFFEFFFFFFFF
     */
    li      a0, -1
    slli    a1, a0, 32
    addi    a2, a1, -1

    /* add mask * FP_MONT_ONE */
    andi    a3, t2, 1
    add     t4, t4, a3
    sltu    t3, t4, a3

    /* reuse ~a1 = 0x00000000FFFFFFFF */
    not     a3, a1
    and     a3, a3, t2
    add     a4, t5, a3
    sltu    a5, a4, t5
    add     t5, a4, t3
    sltu    a6, t5, t3
    or      t3, a5, a6

    add     t0, t0, t3
    sltu    t3, t0, t3

    li      a3, 1
    slli    a3, a3, 32
    and     a3, a3, t2
    add     a4, t1, a3
    sltu    a5, a4, t1
    add     t1, a4, t3
    sltu    a6, t1, t3
    or      t3, a5, a6

    /* =========================================================
     * diagnostic version:
     * keep compensation path, skip ONLY final COND_SUB_P
     * ========================================================= */
    COND_SUB_P t4, t5, t0, t1, a3, a4, a5, a6, t2, t3 

    ld      a7, 32(sp)
    sd      t4,  0(a7)
    sd      t5,  8(a7)
    sd      t0, 16(a7)
    sd      t1, 24(a7)

    ld      s0,  0(sp)
    ld      s1,  8(sp)
    ld      s2, 16(sp)
    ld      s3, 24(sp)
    addi    sp, sp, 40
    ret

.size fp_mont_sqr, .-fp_mont_sqr
    .option pop