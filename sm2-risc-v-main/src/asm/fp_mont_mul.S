    .text
    .align  2
    .globl  fp_mont_mul
    .type   fp_mont_mul, @function

/* ============================================================
 * fp_mont_mul(fp_t* r, const fp_t* a, const fp_t* b)
 *
 * Safe replace version:
 * - 保留原 round0~3 调度
 * - 保留 SM2_REDC_STEP_PRESERVE
 * - 保留双条件减
 * - 只优化 FP_MONT_ONE 补偿段
 * - 条件减中把 not t6, zero 改为 li t6, -1
 * ============================================================ */

    .macro SM2_REDC_STEP_PRESERVE w0, w1, w2, w3, w4, w5
        mv      t6, \w0
        slli    a3, t6, 32
        srli    a4, t6, 32

        mv      \w0, zero

        add     a5, \w1, t6
        sltu    a6, a5, \w1
        sltu    a7, a5, a3
        sub     \w1, a5, a3

        add     a5, \w2, a6
        sltu    a6, a5, \w2
        sltu    a0, a5, a4
        sub     a5, a5, a4
        sltu    a1, a5, a7
        sub     \w2, a5, a7
        or      a7, a0, a1

        add     a5, \w3, a6
        sltu    a6, a5, \w3
        sltu    a0, a5, a3
        sub     a5, a5, a3
        sltu    a1, a5, a7
        sub     \w3, a5, a7
        or      a7, a0, a1

        add     a5, \w4, t6
        sltu    a0, a5, \w4
        add     a5, a5, a6
        sltu    a1, a5, a6
        or      a6, a0, a1
        sltu    a0, a5, a4
        sub     a5, a5, a4
        sltu    a1, a5, a7
        sub     \w4, a5, a7
        or      a7, a0, a1

        add     a5, \w5, a6
        sub     \w5, a5, a7
    .endm

fp_mont_mul:
    addi    sp, sp, -56
    sd      s0,  0(sp)
    sd      s1,  8(sp)
    sd      s2, 16(sp)
    sd      s3, 24(sp)
    sd      s4, 32(sp)
    sd      s5, 40(sp)
    sd      s6, 48(sp)

    mv      s4, a0

    ld      s0,  0(a1)
    ld      s1,  8(a1)
    ld      s2, 16(a1)
    ld      s3, 24(a1)

    /* 常驻模数 limb 常量，尾部条件减复用 */
    li      s5, 0xFFFFFFFF00000000
    li      s6, 0xFFFFFFFEFFFFFFFF

    mv      t0, zero
    mv      t1, zero
    mv      t2, zero
    mv      t3, zero
    mv      t4, zero
    mv      t5, zero

    /* round 0 */
    ld      a1, 0(a2)

    mul     a3, s0, a1
    mulhu   a4, s0, a1
    mul     a5, s1, a1
    mulhu   a6, s1, a1
    mul     a7, s2, a1
    mulhu   a0, s2, a1

    add     t0, t0, a3
    sltu    a3, t0, a3
    add     a4, a4, a3

    mul     t6, s3, a1
    mulhu   a1, s3, a1

    add     t1, t1, a5
    sltu    a3, t1, a5
    add     t1, t1, a4
    sltu    a4, t1, a4
    add     a3, a3, a4
    add     a6, a6, a3

    add     t2, t2, a7
    sltu    a3, t2, a7
    add     t2, t2, a6
    sltu    a4, t2, a6
    add     a3, a3, a4
    add     a0, a0, a3

    add     t3, t3, t6
    sltu    a3, t3, t6
    add     t3, t3, a0
    sltu    a4, t3, a0
    add     a3, a3, a4
    add     a1, a1, a3

    add     t4, t4, a1
    sltu    a3, t4, a1
    add     t5, t5, a3

    SM2_REDC_STEP_PRESERVE t0, t1, t2, t3, t4, t5

    /* round 1 */
    ld      a1, 8(a2)

    mul     a3, s0, a1
    mulhu   a4, s0, a1
    mul     a5, s1, a1
    mulhu   a6, s1, a1
    mul     a7, s2, a1
    mulhu   a0, s2, a1

    add     t1, t1, a3
    sltu    a3, t1, a3
    add     a4, a4, a3

    mul     t6, s3, a1
    mulhu   a1, s3, a1

    add     t2, t2, a5
    sltu    a3, t2, a5
    add     t2, t2, a4
    sltu    a4, t2, a4
    add     a3, a3, a4
    add     a6, a6, a3

    add     t3, t3, a7
    sltu    a3, t3, a7
    add     t3, t3, a6
    sltu    a4, t3, a6
    add     a3, a3, a4
    add     a0, a0, a3

    add     t4, t4, t6
    sltu    a3, t4, t6
    add     t4, t4, a0
    sltu    a4, t4, a0
    add     a3, a3, a4
    add     a1, a1, a3

    add     t5, t5, a1
    sltu    a3, t5, a1
    add     t0, t0, a3

    SM2_REDC_STEP_PRESERVE t1, t2, t3, t4, t5, t0

    /* round 2 */
    ld      a1, 16(a2)

    mul     a3, s0, a1
    mulhu   a4, s0, a1
    mul     a5, s1, a1
    mulhu   a6, s1, a1
    mul     a7, s2, a1
    mulhu   a0, s2, a1

    add     t2, t2, a3
    sltu    a3, t2, a3
    add     a4, a4, a3

    mul     t6, s3, a1
    mulhu   a1, s3, a1

    add     t3, t3, a5
    sltu    a3, t3, a5
    add     t3, t3, a4
    sltu    a4, t3, a4
    add     a3, a3, a4
    add     a6, a6, a3

    add     t4, t4, a7
    sltu    a3, t4, a7
    add     t4, t4, a6
    sltu    a4, t4, a6
    add     a3, a3, a4
    add     a0, a0, a3

    add     t5, t5, t6
    sltu    a3, t5, t6
    add     t5, t5, a0
    sltu    a4, t5, a0
    add     a3, a3, a4
    add     a1, a1, a3

    add     t0, t0, a1
    sltu    a3, t0, a1
    add     t1, t1, a3

    SM2_REDC_STEP_PRESERVE t2, t3, t4, t5, t0, t1

    /* round 3 */
    ld      a1, 24(a2)

    mul     a3, s0, a1
    mulhu   a4, s0, a1
    mul     a5, s1, a1
    mulhu   a6, s1, a1
    mul     a7, s2, a1
    mulhu   a0, s2, a1

    add     t3, t3, a3
    sltu    a3, t3, a3
    add     a4, a4, a3

    mul     t6, s3, a1
    mulhu   a1, s3, a1

    add     t4, t4, a5
    sltu    a3, t4, a5
    add     t4, t4, a4
    sltu    a4, t4, a4
    add     a3, a3, a4
    add     a6, a6, a3

    add     t5, t5, a7
    sltu    a3, t5, a7
    add     t5, t5, a6
    sltu    a4, t5, a6
    add     a3, a3, a4
    add     a0, a0, a3

    add     t0, t0, t6
    sltu    a3, t0, t6
    add     t0, t0, a0
    sltu    a4, t0, a0
    add     a3, a3, a4
    add     a1, a1, a3

    add     t1, t1, a1
    sltu    a3, t1, a1
    add     t2, t2, a3

    SM2_REDC_STEP_PRESERVE t3, t4, t5, t0, t1, t2

    /* 结果：丢弃 t3，r=[t4,t5,t0,t1]，extra=t2 */

    /* ====================================================
     * FP_MONT_ONE 补偿（替换为 sqr 版安全写法）
     * ==================================================== */
    andi    a6, t2, 1           /* bit */
    mv      a5, a6
    neg     t2, a6              /* mask = -bit */

    add     t4, t4, a6
    sltu    t3, t4, a6

    srli    a3, t2, 32          /* mask>>32 = limb1补偿 */
    add     a4, t5, a3
    sltu    a6, a4, t5
    add     t5, a4, t3
    sltu    a7, t5, t3
    or      t3, a6, a7

    add     t0, t0, t3
    sltu    t3, t0, t3

    slli    a3, a5, 32          /* bit<<32 = limb3补偿 */
    add     a4, t1, a3
    sltu    a6, a4, t1
    add     t1, a4, t3
    sltu    a7, t1, t3
    or      t3, a6, a7

    /* ====================================================
     * 条件减 p #1
     * ==================================================== */
    mv      t2, zero
    li      t6, -1

    sub     t3, t4, t6
    sltu    a0, t4, t6
    sub     a3, t3, t2
    sltu    t2, t3, t2
    or      t2, t2, a0

    sub     t3, t5, s5
    sltu    a0, t5, s5
    sub     a4, t3, t2
    sltu    t2, t3, t2
    or      t2, t2, a0

    sub     t3, t0, t6
    sltu    a0, t0, t6
    sub     a5, t3, t2
    sltu    t2, t3, t2
    or      t2, t2, a0

    sub     t3, t1, s6
    sltu    a0, t1, s6
    sub     a6, t3, t2
    sltu    t2, t3, t2
    or      t2, t2, a0

    xori    t2, t2, 1
    sub     t2, x0, t2

    and     a3, a3, t2
    andn    t4, t4, t2
    or      t4, t4, a3

    and     a4, a4, t2
    andn    t5, t5, t2
    or      t5, t5, a4

    and     a5, a5, t2
    andn    t0, t0, t2
    or      t0, t0, a5

    and     a6, a6, t2
    andn    t1, t1, t2
    or      t1, t1, a6

    /* ====================================================
     * 条件减 p #2
     * ==================================================== */
    mv      t2, zero

    sub     t3, t4, t6
    sltu    a0, t4, t6
    sub     a3, t3, t2
    sltu    t2, t3, t2
    or      t2, t2, a0

    sub     t3, t5, s5
    sltu    a0, t5, s5
    sub     a4, t3, t2
    sltu    t2, t3, t2
    or      t2, t2, a0

    sub     t3, t0, t6
    sltu    a0, t0, t6
    sub     a5, t3, t2
    sltu    t2, t3, t2
    or      t2, t2, a0

    sub     t3, t1, s6
    sltu    a0, t1, s6
    sub     a6, t3, t2
    sltu    t2, t3, t2
    or      t2, t2, a0

    xori    t2, t2, 1
    sub     t2, x0, t2

    and     a3, a3, t2
    andn    t4, t4, t2
    or      t4, t4, a3

    and     a4, a4, t2
    andn    t5, t5, t2
    or      t5, t5, a4

    and     a5, a5, t2
    andn    t0, t0, t2
    or      t0, t0, a5

    and     a6, a6, t2
    andn    t1, t1, t2
    or      t1, t1, a6

    sd      t4,  0(s4)
    sd      t5,  8(s4)
    sd      t0, 16(s4)
    sd      t1, 24(s4)

    ld      s0,  0(sp)
    ld      s1,  8(sp)
    ld      s2, 16(sp)
    ld      s3, 24(sp)
    ld      s4, 32(sp)
    ld      s5, 40(sp)
    ld      s6, 48(sp)
    addi    sp, sp, 56
    ret

.size fp_mont_mul, .-fp_mont_mul