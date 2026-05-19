    .text
    .align 4
    .globl fp_mont_inv_pornin_asm
    .type  fp_mont_inv_pornin_asm, @function

#ifdef PORNIN_FULL_INV_PROFILE
    .macro PROFILE_ADD sym, start, tmp0, tmp1, tmp2
    rdcycle \tmp0
    sub     \tmp0, \tmp0, \start
    la      \tmp1, \sym
    ld      \tmp2, 0(\tmp1)
    add     \tmp2, \tmp2, \tmp0
    sd      \tmp2, 0(\tmp1)
    .endm

    .macro PROFILE_INC sym, tmp0, tmp1
    la      \tmp0, \sym
    ld      \tmp1, 0(\tmp0)
    addi    \tmp1, \tmp1, 1
    sd      \tmp1, 0(\tmp0)
    .endm
#endif

/*
 * Self-contained macro set for the fused full inversion path.
 * The split asm files remain available for unit tests and A/B benchmarks,
 * but this file deliberately does not include or call them.
 */

    .macro COND_SUB_SM2_P
    sub     a0, s0, t0
    sltu    a4, s0, t0

    sub     a1, s1, t1
    sltu    a5, s1, t1
    sltu    a6, a1, a4
    sub     a1, a1, a4
    or      a4, a5, a6

    sub     a2, s2, t2
    sltu    a5, s2, t2
    sltu    a6, a2, a4
    sub     a2, a2, a4
    or      a4, a5, a6

    sub     a3, s3, t3
    sltu    a5, s3, t3
    sltu    a6, a3, a4
    sub     a3, a3, a4
    or      a4, a5, a6

    xori    a4, a4, 1
    sub     a4, x0, a4
    xori    t4, a4, -1

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

    .macro LOAD_SM2_P_CONSTS
    li      t0, -1
    li      t1, 0xFFFFFFFF00000000
    li      t2, -1
    li      t3, 0xFFFFFFFEFFFFFFFF
    .endm

/*
 * Fast negation modulo SM2 p for temporary u/v values.
 *
 * If src == 0, this writes p instead of 0.  That representative is equivalent
 * modulo p and has passed the full inversion randomized tests.
 *
 * Requires t0..t3 to contain SM2 p limbs, via LOAD_SM2_P_CONSTS.
 */
    .macro NEG_MOD_SM2_FAST dst, src
    ld      a2, 0(\src)
    ld      a3, 8(\src)
    ld      a4, 16(\src)
    ld      a5, 24(\src)

    sub     t5, t0, a2
    sltu    t6, t0, a2

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

    sd      t5, 0(\dst)
    sd      a6, 8(\dst)
    sd      a7, 16(\dst)
    sd      a0, 24(\dst)
    .endm

/*
 * Register-argument modular lincomb.
 *
 * Input:
 *   s6 = output fp_t*
 *   a1 = c1
 *   a2 = x
 *   a3 = neg_x
 *   a4 = c2
 *   a5 = y
 *   a6 = neg_y
 *
 * Clobbers:
 *   s0-s6, a0-a7, t0-t6
 */
    .macro LINCOMB_REGARGS
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

    ld      t0, 0(s0)
    ld      t1, 0(s1)
    mul     t2, s2, t0
    mulhu   t3, s2, t0
    mul     t4, s3, t1
    mulhu   t5, s3, t1
    ld      a6, 8(s0)
    ld      a7, 8(s1)
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     s4, t3, t5
    add     s4, s4, a0
    mv      a2, t6             /* z0 */

    mul     t2, s2, a6
    mulhu   t3, s2, a6
    mul     t4, s3, a7
    mulhu   t5, s3, a7
    ld      t0, 16(s0)
    ld      t1, 16(s1)
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     t6, t6, s4
    sltu    a1, t6, s4
    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1
    mv      a3, t6             /* z1 */

    mul     t2, s2, t0
    mulhu   t3, s2, t0
    mul     t4, s3, t1
    mulhu   t5, s3, t1
    ld      a6, 24(s0)
    ld      a7, 24(s1)
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     t6, t6, s4
    sltu    a1, t6, s4
    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1
    mv      a4, t6             /* z2 */

    mul     t2, s2, a6
    mulhu   t3, s2, a6
    mul     t4, s3, a7
    mulhu   t5, s3, a7
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     t6, t6, s4
    sltu    a1, t6, s4
    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1
    mv      a5, t6             /* z3, s4 = z4 */

    mv      s0, a2
    mv      s1, a3
    mv      s2, a4
    mv      s3, a5

    add     s0, s0, s4
    sltu    t1, s0, s4
    slli    t2, s4, 32
    srli    t3, s4, 32
    sub     t4, t2, s4
    sltu    t5, t2, s4
    sub     t6, t3, t5
    add     s1, s1, t4
    sltu    a1, s1, t4
    add     s1, s1, t1
    sltu    a2, s1, t1
    add     a3, t6, a1
    add     a3, a3, a2
    add     s2, s2, t3
    sltu    a5, s2, t3
    add     s2, s2, a3
    sltu    a6, s2, a3
    add     a7, a5, a6
    add     s3, s3, t2
    sltu    a1, s3, t2
    add     s3, s3, a7
    sltu    a2, s3, a7
    add     s4, t3, a1
    add     s4, s4, a2

    /*
     * Second fold.  After the first fold, h is only the carry-out from
     * folding a 31-bit-coefficient 320-bit lincomb, so h < 2^32.  This
     * lets us specialize h >> 32 = 0 and ((h << 32) - h) has no borrow.
     */
    add     s0, s0, s4
    sltu    t1, s0, s4
    slli    t2, s4, 32
    sub     t4, t2, s4
    add     s1, s1, t4
    sltu    a1, s1, t4
    add     s1, s1, t1
    sltu    a2, s1, t1
    add     a3, a1, a2
    add     s2, s2, a3
    sltu    a7, s2, a3
    add     s3, s3, t2
    add     s3, s3, a7

    sd      s0, 0(s6)
    sd      s1, 8(s6)
    sd      s2, 16(s6)
    sd      s3, 24(s6)
    .endm

/*
 * Tail-only lincomb entry:
 *   r = c1*u + c2*v mod p
 *
 * Inputs:
 *   s6  = output fp_t*
 *   s7  = c1
 *   s11 = c2
 *   s9  = current slot base, with u at +64 and v at +96
 *   UV_NEG_U / UV_NEG_V already contain -u/-v.
 *
 * This is the same core as LINCOMB_REGARGS, but specialized for the single
 * tail_uv call so it avoids setting up six generic argument registers.
 */
    .macro LINCOMB_TAIL_UV
    srai    t0, s7, 63
    addi    a2, s9, 64
    addi    a3, sp, UV_NEG_U
    xor     t1, a2, a3
    and     t1, t1, t0
    xor     s0, a2, t1
    xor     s2, s7, t0
    sub     s2, s2, t0

    srai    t0, s11, 63
    addi    a5, s9, 96
    addi    a6, sp, UV_NEG_V
    xor     t1, a5, a6
    and     t1, t1, t0
    xor     s1, a5, t1
    xor     s3, s11, t0
    sub     s3, s3, t0

    li      s4, 0

    ld      t0, 0(s0)
    ld      t1, 0(s1)
    mul     t2, s2, t0
    mulhu   t3, s2, t0
    mul     t4, s3, t1
    mulhu   t5, s3, t1
    ld      a6, 8(s0)
    ld      a7, 8(s1)
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     s4, t3, t5
    add     s4, s4, a0
    mv      a2, t6

    mul     t2, s2, a6
    mulhu   t3, s2, a6
    mul     t4, s3, a7
    mulhu   t5, s3, a7
    ld      t0, 16(s0)
    ld      t1, 16(s1)
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     t6, t6, s4
    sltu    a1, t6, s4
    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1
    mv      a3, t6

    mul     t2, s2, t0
    mulhu   t3, s2, t0
    mul     t4, s3, t1
    mulhu   t5, s3, t1
    ld      a6, 24(s0)
    ld      a7, 24(s1)
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     t6, t6, s4
    sltu    a1, t6, s4
    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1
    mv      a4, t6

    mul     t2, s2, a6
    mulhu   t3, s2, a6
    mul     t4, s3, a7
    mulhu   t5, s3, a7
    add     t6, t2, t4
    sltu    a0, t6, t2
    add     t6, t6, s4
    sltu    a1, t6, s4
    add     s4, t3, t5
    add     s4, s4, a0
    add     s4, s4, a1
    mv      a5, t6

    mv      s0, a2
    mv      s1, a3
    mv      s2, a4
    mv      s3, a5

    add     s0, s0, s4
    sltu    t1, s0, s4
    slli    t2, s4, 32
    srli    t3, s4, 32
    sub     t4, t2, s4
    sltu    t5, t2, s4
    sub     t6, t3, t5
    add     s1, s1, t4
    sltu    a1, s1, t4
    add     s1, s1, t1
    sltu    a2, s1, t1
    add     a3, t6, a1
    add     a3, a3, a2
    add     s2, s2, t3
    sltu    a5, s2, t3
    add     s2, s2, a3
    sltu    a6, s2, a3
    add     a7, a5, a6
    add     s3, s3, t2
    sltu    a1, s3, t2
    add     s3, s3, a7
    sltu    a2, s3, a7
    add     s4, t3, a1
    add     s4, s4, a2

    add     s0, s0, s4
    sltu    t1, s0, s4
    slli    t2, s4, 32
    sub     t4, t2, s4
    add     s1, s1, t4
    sltu    a1, s1, t4
    add     s1, s1, t1
    sltu    a2, s1, t1
    add     a3, a1, a2
    add     s2, s2, a3
    sltu    a7, s2, a3
    add     s3, s3, t2
    add     s3, s3, a7

    sd      s0, 0(s6)
    sd      s1, 8(s6)
    sd      s2, 16(s6)
    sd      s3, 24(s6)
    .endm

    .equ UV_NEG_U, 0
    .equ UV_NEG_V, 32
    .equ SCRATCH_T, 64
    .equ F0_SLOT, 208
    .equ G0_SLOT, 216
    .equ F1_SLOT, 224
    .equ G1_SLOT, 232
    .equ SAVE_RA, 256
    .equ SAVE_S7, 264
    .equ SAVE_S8, 272
    .equ SAVE_S9, 280
    .equ SAVE_S10, 288
    .equ SAVE_S11, 296
    .equ SAVE_S0, 304
    .equ SAVE_S1, 312
    .equ SAVE_S2, 320
    .equ SAVE_S3, 328
    .equ SAVE_S4, 336
    .equ SAVE_S5, 344
    .equ SAVE_S6, 352
    .equ SLOT0_A, 384
    .equ SLOT0_B, 416
    .equ SLOT0_U, 448
    .equ SLOT0_V, 480
    .equ SLOT1_A, 512
    .equ SLOT1_B, 544
    .equ SLOT1_U, 576
    .equ SLOT1_V, 608

/*
 * out = bitlen(fp at base), clobbers t0-t4.
 * Uses clz; SpacemiT X60 supports the target CPU tuning used by this project.
 */
    .macro BITLEN_U256 base, out
    ld      t0, 24(\base)
    bnez    t0, 1f
    ld      t0, 16(\base)
    bnez    t0, 2f
    ld      t0, 8(\base)
    bnez    t0, 3f
    ld      t0, 0(\base)
    bnez    t0, 4f
    li      \out, 0
    j       5f
1:
    clz     t1, t0
    li      \out, 256
    sub     \out, \out, t1
    j       5f
2:
    clz     t1, t0
    li      \out, 192
    sub     \out, \out, t1
    j       5f
3:
    clz     t1, t0
    li      \out, 128
    sub     \out, \out, t1
    j       5f
4:
    clz     t1, t0
    li      \out, 64
    sub     \out, \out, t1
5:
    .endm

/*
 * out = low64((fp at base) >> shift). Clobbers t0-t6.
 */
    .macro SHR_LOW64 base, shift, out
    srli    t0, \shift, 6      /* limb */
    andi    t1, \shift, 63     /* off */

    beqz    t0, 1f
    li      t2, 1
    beq     t0, t2, 2f
    li      t2, 2
    beq     t0, t2, 3f
    ld      t3, 24(\base)
    li      t4, 0
    j       4f
1:
    ld      t3, 0(\base)
    ld      t4, 8(\base)
    j       4f
2:
    ld      t3, 8(\base)
    ld      t4, 16(\base)
    j       4f
3:
    ld      t3, 16(\base)
    ld      t4, 24(\base)
4:
    srl     \out, t3, t1
    beqz    t1, 5f
    li      t5, 64
    sub     t5, t5, t1
    sll     t6, t4, t5
    or      \out, \out, t6
5:
    .endm

/*
 * Inline extract_approx64(a,b).
 * Inputs: s9 = current slot base.
 * Outputs: a0 = a_hat, a1 = b_hat.
 */
    .macro EXTRACT_APPROX64_INLINE
    addi    a0, s9, 0
    addi    a1, s9, 32
    BITLEN_U256 a0, a2
    BITLEN_U256 a1, a3
    bge     a2, a3, 1f
    mv      a2, a3
1:
    li      a3, 64
    bge     a2, a3, 2f
    li      a2, 64
2:
    addi    a2, a2, -33       /* shift */

    SHR_LOW64 a0, a2, a4
    li      a5, 0x1FFFFFFFF
    and     a4, a4, a5
    ld      a6, 0(a0)
    li      a7, 0x7FFFFFFF
    and     a6, a6, a7
    slli    a4, a4, 31
    or      a4, a4, a6
    mv      a0, a4

    SHR_LOW64 a1, a2, a4
    li      a5, 0x1FFFFFFFF
    and     a4, a4, a5
    ld      a6, 0(a1)
    li      a7, 0x7FFFFFFF
    and     a6, a6, a7
    slli    a4, a4, 31
    or      a4, a4, a6
    mv      a1, a4
    .endm

/*
 * Inline inner31, inputs a0/a1 and outputs f/g slots.
 */
    .macro INNER31_STEP
    andi    a2, a0, 1
    sltu    a3, a0, a1
    sub     a5, x0, a2
    and     a3, a3, a2
    sub     a4, x0, a3

    xor     t5, a0, a1
    xor     a6, t0, t2
    xor     a7, t1, t3
    and     t5, t5, a4
    and     a6, a6, a4
    and     a7, a7, a4
    xor     a0, a0, t5
    xor     a1, a1, t5
    xor     t0, t0, a6
    xor     t2, t2, a6
    xor     t1, t1, a7
    xor     t3, t3, a7

    srli    a6, a0, 1
    sub     a7, a0, a1
    srli    a7, a7, 1
    xor     t5, a6, a7
    and     t5, t5, a5
    xor     a0, a6, t5

    sub     a6, t0, t2
    sub     a7, t1, t3
    xor     t5, t0, a6
    and     t5, t5, a5
    xor     t0, t0, t5
    xor     t5, t1, a7
    and     t5, t5, a5
    xor     t1, t1, t5

    add     t2, t2, t2
    add     t3, t3, t3
    .endm

    .macro INNER31_INLINE
    li      t0, 1
    li      t1, 0
    li      t2, 0
    li      t3, 1
#ifdef USE_PORNIN_INNER31_HALF_UNROLL
    li      t4, 7
1:
    INNER31_STEP
    INNER31_STEP
    INNER31_STEP
    INNER31_STEP
    addi    t4, t4, -1
    bnez    t4, 1b
    INNER31_STEP
    INNER31_STEP
    INNER31_STEP
#else
    .rept   31
    INNER31_STEP
    .endr
#endif

    mv      a4, t0             /* f0 */
    mv      a5, t1             /* g0 */
    mv      a6, t2             /* f1 */
    mv      a7, t3             /* g1 */
    .endm

    .macro COND_NEG_SLOT slot
    ld      t0, \slot(sp)
    sub     t1, x0, t0
    xor     t2, t0, t1
    and     t2, t2, a7
    xor     t0, t0, t2
    sd      t0, \slot(sp)
    .endm

    .macro COND_NEG_REG reg
    sub     t0, x0, \reg
    xor     t1, \reg, t0
    and     t1, t1, a7
    xor     \reg, \reg, t1
    .endm

/*
 * One limb of two-way interleaved update_ab.
 *
 * Current slot:
 *   s9  = current slot base
 *   s10 = next slot base
 *
 * Coefficients:
 *   a4 = f0
 *   a5 = g0
 *   a6 = f1
 *   a7 = g1
 *
 * Carries:
 *   s4 = carryA
 *   s5 = carryB
 *
 * Output:
 *   zA_i -> zareg
 *   zB_i -> bscratch(sp)
 *
 * This computes one limb of:
 *   zA = f0*a + g0*b
 *   zB = f1*a + g1*b
 */
    .macro AB2_LIMB_FULL aoff, boff, zareg, bscratch
    /*
     * Load common input limbs.
     * a limb from current a, b limb from current b.
     */
    ld      t0, \aoff(s9)       /* ai */
    ld      t1, \boff(s9)       /* bi */

    /*
     * Start A and B multiplications using ai.
     * A lane: f0 * ai
     * B lane: f1 * ai
     */
    mul     t2, a4, t0          /* A lo1 */
    mulhsu  t3, a4, t0          /* A hi1, signed*unsigned */
    mul     t4, a6, t0          /* B lo1 */
    mulhsu  t5, a6, t0          /* B hi1, signed*unsigned */

    /*
     * Start A and B multiplications using bi.
     * A lane: g0 * bi
     * B lane: g1 * bi
     */
    mul     t6, a5, t1          /* A lo2 */
    mulhsu  a0, a5, t1          /* A hi2 */
    mul     a1, a7, t1          /* B lo2 */
    mulhsu  t0, a7, t1          /* B hi2 */

    /*
     * Finish A low limb and signed carry.
     */
    add     a2, t2, t6          /* A low = lo1 + lo2 */
    sltu    a3, a2, t2          /* carry0 */
    add     a2, a2, s4          /* + old carryA low */
    sltu    t6, a2, s4          /* carry1 */
    srai    t1, s4, 63          /* sign_extend(old carryA) */

    add     s4, t3, a0
    add     s4, s4, a3
    add     s4, s4, t6
    add     s4, s4, t1          /* next carryA */

    /*
     * Finish B low limb and signed carry.
     */
    add     a0, t4, a1          /* B low = lo1 + lo2 */
    sltu    t4, a0, t4          /* carry0 */
    add     a0, a0, s5          /* + old carryB low */
    sltu    a1, a0, s5          /* carry1 */
    srai    t1, s5, 63          /* sign_extend(old carryB) */

    add     s5, t5, t0
    add     s5, s5, t4
    add     s5, s5, a1
    add     s5, s5, t1          /* next carryB */

    mv      \zareg, a2          /* zA_i */
    sd      a0, \bscratch(sp)   /* zB_i */
    .endm

    .macro AB2_LIMB_FIRST_PREFETCH aoff, boff, nextaoff, nextboff, zareg, bscratch
    ld      t0, \aoff(s9)
    ld      t1, \boff(s9)

    mul     t2, a4, t0
    mulhsu  t3, a4, t0
    mul     t4, a6, t0
    mulhsu  t5, a6, t0

    mul     t6, a5, t1
    mulhsu  a0, a5, t1
    mul     a1, a7, t1
    mulhsu  t0, a7, t1

    add     a2, t2, t6
    sltu    a3, a2, t2
    add     s4, t3, a0
    add     s4, s4, a3

    ld      s6, \nextaoff(s9)
    ld      a3, \nextboff(s9)

    add     a0, t4, a1
    sltu    t4, a0, t4
    add     s5, t5, t0
    add     s5, s5, t4

    mv      \zareg, a2
    sd      a0, \bscratch(sp)
    .endm

    .macro AB2_LIMB_FULL_PREFETCH aoff, boff, nextaoff, nextboff, zareg, bscratch
    ld      t0, \aoff(s9)
    ld      t1, \boff(s9)

    mul     t2, a4, t0
    mulhsu  t3, a4, t0
    mul     t4, a6, t0
    mulhsu  t5, a6, t0

    mul     t6, a5, t1
    mulhsu  a0, a5, t1
    mul     a1, a7, t1
    mulhsu  t0, a7, t1

    add     a2, t2, t6
    sltu    a3, a2, t2
    add     a2, a2, s4
    sltu    t6, a2, s4
    srai    t1, s4, 63

    add     s4, t3, a0
    add     s4, s4, a3
    ld      s6, \nextaoff(s9)
    ld      a3, \nextboff(s9)
    add     s4, s4, t6
    add     s4, s4, t1

    add     a0, t4, a1
    sltu    t4, a0, t4
    add     a0, a0, s5
    sltu    a1, a0, s5
    srai    t1, s5, 63

    add     s5, t5, t0
    add     s5, s5, t4
    add     s5, s5, a1
    add     s5, s5, t1

    mv      \zareg, a2
    sd      a0, \bscratch(sp)
    .endm

    .macro AB2_LIMB_LOADED_PREFETCH nextaoff, nextboff, zareg, bscratch
    mul     t2, a4, s6
    mulhsu  t3, a4, s6
    mul     t4, a6, s6
    mulhsu  t5, a6, s6

    mul     t6, a5, a3
    mulhsu  a0, a5, a3
    mul     a1, a7, a3
    mulhsu  t0, a7, a3

    add     a2, t2, t6
    sltu    a3, a2, t2
    add     a2, a2, s4
    sltu    t6, a2, s4
    srai    t1, s4, 63

    add     s4, t3, a0
    add     s4, s4, a3
    ld      s6, \nextaoff(s9)
    ld      a3, \nextboff(s9)
    add     s4, s4, t6
    add     s4, s4, t1

    add     a0, t4, a1
    sltu    t4, a0, t4
    add     a0, a0, s5
    sltu    a1, a0, s5
    srai    t1, s5, 63

    add     s5, t5, t0
    add     s5, s5, t4
    add     s5, s5, a1
    add     s5, s5, t1

    mv      \zareg, a2
    sd      a0, \bscratch(sp)
    .endm

    .macro AB2_LIMB_LOADED_LAST zareg, bscratch
    mul     t2, a4, s6
    mulhsu  t3, a4, s6
    mul     t4, a6, s6
    mulhsu  t5, a6, s6

    mul     t6, a5, a3
    mulhsu  a0, a5, a3
    mul     a1, a7, a3
    mulhsu  t0, a7, a3

    add     a2, t2, t6
    sltu    a3, a2, t2
    add     a2, a2, s4
    sltu    t6, a2, s4
    srai    t1, s4, 63

    add     s4, t3, a0
    add     s4, s4, a3
    add     s4, s4, t6
    add     s4, s4, t1

    add     a0, t4, a1
    sltu    t4, a0, t4
    add     a0, a0, s5
    sltu    a1, a0, s5
    srai    t1, s5, 63

    add     s5, t5, t0
    add     s5, s5, t4
    add     s5, s5, a1
    add     s5, s5, t1

    mv      \zareg, a2
    sd      a0, \bscratch(sp)
    .endm

/*
 * abs(z0..z4) >> 31 and store result.
 *
 * Input:
 *   z0,z1,z2,z3,z4 are registers.
 *   outreg is pointer to output fp_t.
 *
 * Output:
 *   a7 = mask, where mask = -1 if original z was negative, else 0.
 *
 * Clobbers:
 *   t0-t4, a0-a3, a6, a7
 */
    .macro AB_ABS_SHIFT_STORE_MASK z0, z1, z2, z3, z4, outreg
    srli    a7, \z4, 63        /* neg flag */
    sub     a6, x0, a7         /* mask = -neg */
    mv      a0, a7             /* add_one = neg */

    xor     t0, \z0, a6
    add     \z0, t0, a0
    sltu    a0, \z0, t0

    xor     t0, \z1, a6
    add     \z1, t0, a0
    sltu    a0, \z1, t0

    xor     t0, \z2, a6
    add     \z2, t0, a0
    sltu    a0, \z2, t0

    xor     t0, \z3, a6
    add     \z3, t0, a0
    sltu    a0, \z3, t0

    xor     t0, \z4, a6
    add     \z4, t0, a0

    /*
     * out[0] = (z0 >> 31) | (z1 << 33)
     */
    srli    t1, \z0, 31
    slli    t2, \z1, 33
    or      t1, t1, t2
    sd      t1, 0(\outreg)

    srli    t1, \z1, 31
    slli    t2, \z2, 33
    or      t1, t1, t2
    sd      t1, 8(\outreg)

    srli    t1, \z2, 31
    slli    t2, \z3, 33
    or      t1, t1, t2
    sd      t1, 16(\outreg)

    srli    t1, \z3, 31
    slli    t2, \z4, 33
    or      t1, t1, t2
    sd      t1, 24(\outreg)

    /*
     * Leave a7 as mask for COND_NEG_SLOT.
     */
    mv      a7, a6
    .endm


/*
 * Inline update_ab, two-way interleaved version.
 *
 * Current slot base:
 *   s9
 *
 * Next slot base:
 *   s10
 *
 * Output:
 *   new_a -> s10 + 0
 *   new_b -> s10 + 32
 *
 * This replaces the old sequential per-output lincomb with a two-way
 * interleaved update for new_a/new_b.
 */
    .macro UPDATE_AB_INLINE
    /*
     * Clear signed carries:
     *   s4 = carryA
     *   s5 = carryB
     */
    li      s4, 0
    li      s5, 0

    /*
     * Accumulate both lanes.
     *
     * A z limbs:
     *   s0,s1,s2,s3,s4
     *
     * B z limbs:
     *   SCRATCH_T + 0,8,16,24 and final s5
     *
     * current slot layout:
     *   a at s9 + 0
     *   b at s9 + 32
     */
    AB2_LIMB_FIRST_PREFETCH 0, 32, 8, 40, s0, SCRATCH_T+0
    AB2_LIMB_LOADED_PREFETCH 16, 48, s1, SCRATCH_T+8
    AB2_LIMB_LOADED_PREFETCH 24, 56, s2, SCRATCH_T+16
    AB2_LIMB_LOADED_LAST s3, SCRATCH_T+24
    sd      s5, SCRATCH_T+32(sp)   /* zB4 */
    sd      a6, F1_SLOT(sp)        /* preserve f1 across A abs */
    sd      a7, G1_SLOT(sp)        /* preserve g1 across A abs */

    /*
     * Store new_a and conditionally negate f0/g0.
     */
    addi    s6, s10, 0
    AB_ABS_SHIFT_STORE_MASK s0, s1, s2, s3, s4, s6

    /*
     * Start loading B lane before fixing f0/g0.  These loads are independent
     * and can overlap with the conditional coefficient negation on in-order
     * dual-issue cores.
     */
    ld      s0, SCRATCH_T+0(sp)
    ld      s1, SCRATCH_T+8(sp)
    ld      s2, SCRATCH_T+16(sp)
    ld      s3, SCRATCH_T+24(sp)
    ld      s4, SCRATCH_T+32(sp)

    COND_NEG_REG a4
    COND_NEG_REG a5
    sd      a4, F0_SLOT(sp)
    sd      a5, G0_SLOT(sp)

    /*
     * Store new_b, then conditionally negate f1/g1.
     */
    addi    s6, s10, 32
    AB_ABS_SHIFT_STORE_MASK s0, s1, s2, s3, s4, s6
    ld      a4, F1_SLOT(sp)
    ld      a5, G1_SLOT(sp)
    COND_NEG_REG a4
    COND_NEG_REG a5
    sd      a4, F1_SLOT(sp)
    sd      a5, G1_SLOT(sp)
    .endm
/*
 * Inline update_uv:
 *   current slot base = s9, next slot base = s10
 */
    .macro UPDATE_UV_INLINE
    /*
     * neg_u = -u mod p
     * neg_v = -v mod p
     */
    LOAD_SM2_P_CONSTS
    addi    s0, sp, UV_NEG_U
    addi    s1, s9, 64
    NEG_MOD_SM2_FAST s0, s1

    addi    s0, sp, UV_NEG_V
    addi    s1, s9, 96
    NEG_MOD_SM2_FAST s0, s1

    /*
     * new_u = f0*u + g0*v mod p
     */
    addi    s6, s10, 64
    ld      a1, F0_SLOT(sp)
    addi    a2, s9, 64
    addi    a3, sp, UV_NEG_U
    ld      a4, G0_SLOT(sp)
    addi    a5, s9, 96
    addi    a6, sp, UV_NEG_V
    LINCOMB_REGARGS

    /*
     * new_v = f1*u + g1*v mod p
     */
    addi    s6, s10, 96
    ld      a1, F1_SLOT(sp)
    addi    a2, s9, 64
    addi    a3, sp, UV_NEG_U
    ld      a4, G1_SLOT(sp)
    addi    a5, s9, 96
    addi    a6, sp, UV_NEG_V
    LINCOMB_REGARGS
    .endm

    .macro TAIL14_INLINE
    ld      a0, 0(s9)
    ld      a1, 32(s9)
    li      t0, 1
    li      t1, 0
    li      t2, 0
    li      t3, 1
    li      t4, 14
1:
    andi    a2, a0, 1
    sltu    a3, a0, a1
    sub     a5, x0, a2
    and     a3, a3, a2
    sub     a4, x0, a3

    xor     t5, a0, a1
    xor     a6, t0, t2
    xor     a7, t1, t3
    and     t5, t5, a4
    and     a6, a6, a4
    and     a7, a7, a4
    xor     a0, a0, t5
    xor     a1, a1, t5
    xor     t0, t0, a6
    xor     t2, t2, a6
    xor     t1, t1, a7
    xor     t3, t3, a7

    srli    a6, a0, 1
    sub     a7, a0, a1
    srli    a7, a7, 1
    xor     t5, a6, a7
    and     t5, t5, a5
    xor     a0, a6, t5

    sub     a6, t0, t2
    sub     a7, t1, t3
    xor     t5, t0, a6
    and     t5, t5, a5
    xor     t0, t0, t5
    xor     t5, t1, a7
    and     t5, t5, a5
    xor     t1, t1, t5

    add     t2, t2, t2
    add     t3, t3, t3
    addi    t4, t4, -1
    bnez    t4, 1b
    mv      s7, t2             /* F1 */
    mv      s11, t3            /* G1 */
    .endm

/*
 * void fp_mont_inv_pornin_asm(fp_t *r, const fp_t *a_bar);
 *
 * Fully fused experimental Montgomery-domain Pornin inversion:
 *   - inline extract_approx64
 *   - inline inner31
 *   - inline update_ab
 *   - inline update_uv
 *   - inline tail14
 *   - tail_uv uses inline lincomb
 *   - final multiply calls fp_mont_mul
 */
fp_mont_inv_pornin_asm:
    addi    sp, sp, -640
    sd      ra, SAVE_RA(sp)
    sd      s7, SAVE_S7(sp)
    sd      s8, SAVE_S8(sp)
    sd      s9, SAVE_S9(sp)
    sd      s10, SAVE_S10(sp)
    sd      s11, SAVE_S11(sp)
    sd      s0, SAVE_S0(sp)
    sd      s1, SAVE_S1(sp)
    sd      s2, SAVE_S2(sp)
    sd      s3, SAVE_S3(sp)
    sd      s4, SAVE_S4(sp)
    sd      s5, SAVE_S5(sp)
    sd      s6, SAVE_S6(sp)

    mv      s8, a0             /* output r */

    ld      t0, 0(a1)
    ld      t1, 8(a1)
    ld      t2, 16(a1)
    ld      t3, 24(a1)
    or      t4, t0, t1
    or      t4, t4, t2
    or      t4, t4, t3
    bnez    t4, .L_full_nonzero
    sd      x0, 0(s8)
    sd      x0, 8(s8)
    sd      x0, 16(s8)
    sd      x0, 24(s8)
    j       .L_full_done

.L_full_nonzero:
#ifdef PORNIN_FULL_INV_PROFILE
    PROFILE_INC pornin_full_prof_count, t4, t5
#endif

    sd      t0, SLOT0_A(sp)
    sd      t1, SLOT0_A+8(sp)
    sd      t2, SLOT0_A+16(sp)
    sd      t3, SLOT0_A+24(sp)

    la      t4, FP_P
    ld      t0, 0(t4)
    ld      t1, 8(t4)
    ld      t2, 16(t4)
    ld      t3, 24(t4)
    sd      t0, SLOT0_B(sp)
    sd      t1, SLOT0_B+8(sp)
    sd      t2, SLOT0_B+16(sp)
    sd      t3, SLOT0_B+24(sp)

    la      t4, .LPORNIN_FP_R2
    ld      t0, 0(t4)
    ld      t1, 8(t4)
    ld      t2, 16(t4)
    ld      t3, 24(t4)
    sd      t0, SLOT0_U(sp)
    sd      t1, SLOT0_U+8(sp)
    sd      t2, SLOT0_U+16(sp)
    sd      t3, SLOT0_U+24(sp)

    sd      x0, SLOT0_V(sp)
    sd      x0, SLOT0_V+8(sp)
    sd      x0, SLOT0_V+16(sp)
    sd      x0, SLOT0_V+24(sp)

    addi    s9, sp, SLOT0_A
    addi    s10, sp, SLOT1_A
    li      s11, 16

    .align 4
.L_full_round:
#ifdef PORNIN_FULL_INV_PROFILE
    rdcycle s7
#endif
    EXTRACT_APPROX64_INLINE
#ifdef PORNIN_FULL_INV_PROFILE
    PROFILE_ADD pornin_full_prof_extract, s7, t4, t5, t6
    rdcycle s7
#endif
    INNER31_INLINE
#ifdef PORNIN_FULL_INV_PROFILE
    PROFILE_ADD pornin_full_prof_inner31, s7, t4, t5, t6
    rdcycle s7
#endif
    UPDATE_AB_INLINE
#ifdef PORNIN_FULL_INV_PROFILE
    PROFILE_ADD pornin_full_prof_update_ab, s7, t4, t5, t6
    rdcycle s7
#endif
    UPDATE_UV_INLINE
#ifdef PORNIN_FULL_INV_PROFILE
    PROFILE_ADD pornin_full_prof_update_uv, s7, t4, t5, t6
#endif

    mv      s7, s9
    mv      s9, s10
    mv      s10, s7
    addi    s11, s11, -1
    bnez    s11, .L_full_round

    /*
     * After 16 swaps, s9 points back to slot0 (the final current slot).
     */
#ifdef PORNIN_FULL_INV_PROFILE
    rdcycle s10
#endif
    TAIL14_INLINE
#ifdef PORNIN_FULL_INV_PROFILE
    PROFILE_ADD pornin_full_prof_tail14, s10, t4, t5, t6
#endif

    /*
     * tail_uv: next_v = F1*u + G1*v mod p, reusing slot1_v as output.
     */
    LOAD_SM2_P_CONSTS
    addi    s0, sp, UV_NEG_U
    addi    s1, s9, 64
    NEG_MOD_SM2_FAST s0, s1
    addi    s0, sp, UV_NEG_V
    addi    s1, s9, 96
    NEG_MOD_SM2_FAST s0, s1

    /*
     * tail_uv:
     *   next_v = F1*u + G1*v mod p
     */
#ifdef PORNIN_FULL_INV_PROFILE
    rdcycle s10
#endif
    addi    s6, sp, SLOT1_V
    LINCOMB_TAIL_UV
#ifdef PORNIN_FULL_INV_PROFILE
    PROFILE_ADD pornin_full_prof_tail_uv, s10, t4, t5, t6
#endif

    /*
     * final: r = v * C510 mod p.
     */
#ifdef PORNIN_FULL_INV_PROFILE
    rdcycle s7
#endif
    mv      a0, s8
    addi    a1, sp, SLOT1_V
    la      a2, .LPORNIN_FP_INV2_POW_510_BAR
    call    fp_mont_mul
#ifdef PORNIN_FULL_INV_PROFILE
    PROFILE_ADD pornin_full_prof_final_mul, s7, t4, t5, t6
#endif

.L_full_done:
    ld      ra, SAVE_RA(sp)
    ld      s7, SAVE_S7(sp)
    ld      s8, SAVE_S8(sp)
    ld      s9, SAVE_S9(sp)
    ld      s10, SAVE_S10(sp)
    ld      s11, SAVE_S11(sp)
    ld      s0, SAVE_S0(sp)
    ld      s1, SAVE_S1(sp)
    ld      s2, SAVE_S2(sp)
    ld      s3, SAVE_S3(sp)
    ld      s4, SAVE_S4(sp)
    ld      s5, SAVE_S5(sp)
    ld      s6, SAVE_S6(sp)
    addi    sp, sp, 640
    ret

    .size fp_mont_inv_pornin_asm, .-fp_mont_inv_pornin_asm

    .section .rodata
    .align 3
.LPORNIN_FP_R2:
    .dword 0x0000000200000003
    .dword 0x00000002FFFFFFFF
    .dword 0x0000000100000001
    .dword 0x0000000400000002

.LPORNIN_FP_INV2_POW_510_BAR:
    .dword 0xFFFFFFE400000013
    .dword 0xFFFFFFF700000018
    .dword 0xFFFFFFF00000000B
    .dword 0xFFFFFFEF00000017

#ifdef PORNIN_FULL_INV_PROFILE
    .section .bss
    .align 3
    .globl pornin_full_prof_count
    .globl pornin_full_prof_extract
    .globl pornin_full_prof_inner31
    .globl pornin_full_prof_update_ab
    .globl pornin_full_prof_update_uv
    .globl pornin_full_prof_tail14
    .globl pornin_full_prof_tail_uv
    .globl pornin_full_prof_final_mul
pornin_full_prof_count:
    .zero 8
pornin_full_prof_extract:
    .zero 8
pornin_full_prof_inner31:
    .zero 8
pornin_full_prof_update_ab:
    .zero 8
pornin_full_prof_update_uv:
    .zero 8
pornin_full_prof_tail14:
    .zero 8
pornin_full_prof_tail_uv:
    .zero 8
pornin_full_prof_final_mul:
    .zero 8
#endif
