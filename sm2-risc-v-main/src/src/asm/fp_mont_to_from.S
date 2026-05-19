    .text
    .align  2

    .globl  fp_to_mont
    .type   fp_to_mont, @function

    .globl  fp_from_mont
    .type   fp_from_mont, @function

/* ============================================================
 * 复用 fast fp_mont_mul 的 SM2 专用 REDC 步
 * ============================================================ */
        .macro SM2_REDC_STEP_PRESERVE w0, w1, w2, w3, w4, w5
        /* m = current low limb = w0 */
        mv      t6, \w0
        slli    a3, t6, 32
        srli    a4, t6, 32

        /*
         * 这里不能直接把 w0 清零。
         * 当前宏采用“参数轮转”来模拟窗口左移：
         *   REDC(t0,t1,t2,t3,t4,t5)
         *   REDC(t1,t2,t3,t4,t5,t0)
         * 所以本轮最终产生的最高位，必须回填到旧 w0，
         * 供下一轮作为新的最高 limb 使用。
         */

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

        /*
         * limb5 更新后，最终溢出/借位必须写回 w0，
         * 这样下一轮轮转参数时，旧 w0 才能作为新的最高位。
         */
        add     a5, \w5, a6
        sltu    a0, a5, a6
        sub     \w5, a5, a7
        sltu    a1, a5, a7
        or      \w0, a0, a1
    .endm

    .macro MULC2_SAFE lo, hi, x, tmp1, tmp2
        slli    \lo, \x, 32
        srli    \hi, \x, 32
        add     \tmp1, \lo, \x
        sltu    \tmp2, \tmp1, \lo
        mv      \lo, \tmp1
        add     \hi, \hi, \tmp2
    .endm

    .macro MULC3_SAFE lo, hi, x, tmp1, tmp2, tmp3
        slli    \lo, \x, 34
        srli    \hi, \x, 30

        srli    \tmp2, \x, 63
        sh1add  \tmp3, \x, \lo
        sltu    \tmp1, \tmp3, \lo
        mv      \lo, \tmp3

        add     \hi, \hi, \tmp2
        add     \hi, \hi, \tmp1
    .endm


/* ============================================================
 * fp_to_mont(r, a)
 *
 * a -> aR mod p = REDC(a * R^2)
 *
 * 混合特化策略：
 *   c0 = 0x0000000200000003  -> 通用 mul/mulhu
 *   c1 = 0x00000002FFFFFFFF  -> 通用 mul/mulhu
 *   c2 = 0x0000000100000001  -> 特化为 shift+add
 *   c3 = 0x0000000400000002  -> 特化为 shift+add
 * ============================================================ */
fp_to_mont:
    addi    sp, sp, -56
    sd      s0,  0(sp)
    sd      s1,  8(sp)
    sd      s2, 16(sp)
    sd      s3, 24(sp)
    sd      s4, 32(sp)
    sd      s5, 40(sp)
    sd      s6, 48(sp)

    mv      s4, a0

    /* a[0..3] */
    ld      s0,  0(a1)
    ld      s1,  8(a1)
    ld      s2, 16(a1)
    ld      s3, 24(a1)

    /* 常驻模数 limb 常量，供尾部两次条件减法复用 */
    li      s5, 0xFFFFFFFF00000000
    li      s6, 0xFFFFFFFEFFFFFFFF

    mv      t0, zero
    mv      t1, zero
    mv      t2, zero
    mv      t3, zero
    mv      t4, zero
    mv      t5, zero

    /* ====================================================
     * round 0 : window [t0,t1,t2,t3,t4,t5], bi = c0
     * c0 = 0x0000000200000003
     * ==================================================== */
    li      a1, 0x0000000200000003

    mul     a3, s0, a1
    mulhu   a4, s0, a1
    mul     t6, s1, a1
    mulhu   a0, s1, a1

    add     t0, t0, a3
    sltu    a6, t0, a3
    add     a6, a4, a6

    mul     a3, s2, a1
    mulhu   a4, s2, a1

    add     t1, t1, t6
    sltu    a7, t1, t6
    add     t1, t1, a6
    sltu    t6, t1, a6
    add     a7, a7, t6
    add     a6, a0, a7

    mul     t6, s3, a1
    mulhu   a0, s3, a1

    add     t2, t2, a3
    sltu    a7, t2, a3
    add     t2, t2, a6
    sltu    a3, t2, a6
    add     a7, a7, a3
    add     a6, a4, a7

    add     t3, t3, t6
    sltu    a7, t3, t6
    add     t3, t3, a6
    sltu    a3, t3, a6
    add     a7, a7, a3
    add     a6, a0, a7

    add     t4, t4, a6
    sltu    a7, t4, a6
    add     t5, t5, a7

    SM2_REDC_STEP_PRESERVE t0, t1, t2, t3, t4, t5

    /* ====================================================
     * round 1 : window [t1,t2,t3,t4,t5,t0], bi = c1
     * c1 = 0x00000002FFFFFFFF
     * ==================================================== */
    li      a1, 0x00000002FFFFFFFF

    mul     a3, s0, a1
    mulhu   a4, s0, a1
    mul     t6, s1, a1
    mulhu   a0, s1, a1

    add     t1, t1, a3
    sltu    a6, t1, a3
    add     a6, a4, a6

    mul     a3, s2, a1
    mulhu   a4, s2, a1

    add     t2, t2, t6
    sltu    a7, t2, t6
    add     t2, t2, a6
    sltu    t6, t2, a6
    add     a7, a7, t6
    add     a6, a0, a7

    mul     t6, s3, a1
    mulhu   a0, s3, a1

    add     t3, t3, a3
    sltu    a7, t3, a3
    add     t3, t3, a6
    sltu    a3, t3, a6
    add     a7, a7, a3
    add     a6, a4, a7

    add     t4, t4, t6
    sltu    a7, t4, t6
    add     t4, t4, a6
    sltu    a3, t4, a6
    add     a7, a7, a3
    add     a6, a0, a7

    add     t5, t5, a6
    sltu    a7, t5, a6
    add     t0, t0, a7

    SM2_REDC_STEP_PRESERVE t1, t2, t3, t4, t5, t0

    /* ====================================================
     * round 2 : window [t2,t3,t4,t5,t0,t1], bi = c2
     * c2 = 0x0000000100000001
     * ==================================================== */
    MULC2_SAFE a3, a4, s0, a1, a5
    MULC2_SAFE t6, a0, s1, a1, a5

    add     t2, t2, a3
    sltu    a6, t2, a3
    add     a6, a4, a6

    MULC2_SAFE a3, a4, s2, a1, a5

    add     t3, t3, t6
    sltu    a7, t3, t6
    add     t3, t3, a6
    sltu    t6, t3, a6
    add     a7, a7, t6
    add     a6, a0, a7

    MULC2_SAFE t6, a0, s3, a1, a5

    add     t4, t4, a3
    sltu    a7, t4, a3
    add     t4, t4, a6
    sltu    a3, t4, a6
    add     a7, a7, a3
    add     a6, a4, a7

    add     t5, t5, t6
    sltu    a7, t5, t6
    add     t5, t5, a6
    sltu    a3, t5, a6
    add     a7, a7, a3
    add     a6, a0, a7

    add     t0, t0, a6
    sltu    a7, t0, a6
    add     t1, t1, a7

    SM2_REDC_STEP_PRESERVE t2, t3, t4, t5, t0, t1

    /* ====================================================
     * round 3 : window [t3,t4,t5,t0,t1,t2], bi = c3
     * c3 = 0x0000000400000002
     * ==================================================== */
    MULC3_SAFE a3, a4, s0, a1, a5, a2
    MULC3_SAFE t6, a0, s1, a1, a5, a2

    add     t3, t3, a3
    sltu    a6, t3, a3
    add     a6, a4, a6

    MULC3_SAFE a3, a4, s2, a1, a5, a2

    add     t4, t4, t6
    sltu    a7, t4, t6
    add     t4, t4, a6
    sltu    t6, t4, a6
    add     a7, a7, t6
    add     a6, a0, a7

    MULC3_SAFE t6, a0, s3, a1, a5, a2

    add     t5, t5, a3
    sltu    a7, t5, a3
    add     t5, t5, a6
    sltu    a3, t5, a6
    add     a7, a7, a3
    add     a6, a4, a7

    add     t0, t0, t6
    sltu    a7, t0, t6
    add     t0, t0, a6
    sltu    a3, t0, a6
    add     a7, a7, a3
    add     a6, a0, a7

    add     t1, t1, a6
    sltu    a7, t1, a6
    add     t2, t2, a7

    SM2_REDC_STEP_PRESERVE t3, t4, t5, t0, t1, t2

    /* 结果：丢弃 t3，r=[t4,t5,t0,t1]，extra=t2 */

    /* ====================================================
     * FP_MONT_ONE 补偿
     * ==================================================== */
    andi    t2, t2, 1
    sub     t2, x0, t2

    addi    a3, zero, 1
    and     a3, a3, t2
    add     t4, t4, a3
    sltu    t3, t4, a3

    /* 0x00000000FFFFFFFF = ~0xFFFFFFFF00000000 = ~s5 */
    not     a3, s5
    and     a3, a3, t2
    add     a4, t5, a3
    sltu    a5, a4, t5
    add     t5, a4, t3
    sltu    a6, t5, t3
    or      t3, a5, a6

    add     t0, t0, t3
    sltu    t3, t0, t3

    /* 0x0000000100000000 = (s5 >> 32) + 1 */
    srli    a3, s5, 32
    addi    a3, a3, 1
    and     a3, a3, t2
    add     a4, t1, a3
    sltu    a5, a4, t1
    add     t1, a4, t3
    sltu    a6, t1, t3
    or      t3, a5, a6

    /* ====================================================
     * 条件减 p #1
     * ==================================================== */
    mv      t2, zero
    not     t6, zero

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

.size fp_to_mont, .-fp_to_mont


/* ============================================================
 * fp_from_mont(r, a)
 *
 * aR -> a = REDC(aR * 1)
 *
 * 初始窗口：
 *   [a0, a1, a2, a3, 0, 0]
 * 连做 4 轮 REDC
 * ============================================================ */
fp_from_mont:
    addi    sp, sp, -32
    sd      s4,  0(sp)
    sd      s5,  8(sp)
    sd      s6, 16(sp)

    mv      s4, a0

    /* 常驻常量 */
    li      s5, 0xFFFFFFFF00000000
    li      s6, 0xFFFFFFFEFFFFFFFF

    /* 初始窗口 [t0,t1,t2,t3,t4,t5] = [a0,a1,a2,a3,0,0] */
    ld      t0,  0(a1)
    ld      t1,  8(a1)
    ld      t2, 16(a1)
    ld      t3, 24(a1)
    mv      t4, zero
    mv      t5, zero

    /* round 0 : [t0,t1,t2,t3,t4,t5] */
    SM2_REDC_STEP_PRESERVE t0, t1, t2, t3, t4, t5

    /* round 1 : [t1,t2,t3,t4,t5,t0] */
    SM2_REDC_STEP_PRESERVE t1, t2, t3, t4, t5, t0

    /* round 2 : [t2,t3,t4,t5,t0,t1] */
    SM2_REDC_STEP_PRESERVE t2, t3, t4, t5, t0, t1

    /* round 3 : [t3,t4,t5,t0,t1,t2] */
    SM2_REDC_STEP_PRESERVE t3, t4, t5, t0, t1, t2

    /* 结果：丢弃 t3，r=[t4,t5,t0,t1]，extra=t2 */

    /* FP_MONT_ONE 补偿 */
    andi    t2, t2, 1
    sub     t2, x0, t2

    addi    a3, zero, 1
    and     a3, a3, t2
    add     t4, t4, a3
    sltu    t3, t4, a3

    /* 0x00000000FFFFFFFF = ~0xFFFFFFFF00000000 = ~s5 */
    not     a3, s5
    and     a3, a3, t2
    add     a4, t5, a3
    sltu    a5, a4, t5
    add     t5, a4, t3
    sltu    a6, t5, t3
    or      t3, a5, a6

    add     t0, t0, t3
    sltu    t3, t0, t3

    /* 0x0000000100000000 = (s5 >> 32) + 1 */
    srli    a3, s5, 32
    addi    a3, a3, 1
    and     a3, a3, t2
    add     a4, t1, a3
    sltu    a5, a4, t1
    add     t1, a4, t3
    sltu    a6, t1, t3
    or      t3, a5, a6

    /* 条件减 p #1 */
    mv      t2, zero
    not     t6, zero

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

    /* 条件减 p #2 */
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

    ld      s4,  0(sp)
    ld      s5,  8(sp)
    ld      s6, 16(sp)
    addi    sp, sp, 32
    ret

.size fp_from_mont, .-fp_from_mont
