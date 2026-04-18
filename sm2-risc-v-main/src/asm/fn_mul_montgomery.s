# =============================================================================
# fn_mul_montgomery_sched.S
#
# RISC-V RV64 + M extension
#
# Implements:
#   void fn_mul_4x4_u512_ct(uint64_t out[8], const fn_t* a, const fn_t* b);
#
# Notes:
# - out[0..7] = a(4 limbs) * b(4 limbs), little-endian limbs
# - Constant-time style (no data-dependent branches)
# - Uses Comba column accumulation with a 192-bit accumulator
# - Scheduled version: overlaps mul/mulhu issue with ACC chain where possible
# =============================================================================

    .option norvc
    .text
    .align 2

    .macro ACC192_ADD128 acc0, acc1, acc2, lo, hi, tmp1, tmp2
        add     \acc0, \acc0, \lo
        sltu    \tmp1, \acc0, \lo

        add     \acc1, \acc1, \hi
        sltu    \tmp2, \acc1, \hi

        add     \acc1, \acc1, \tmp1
        sltu    \tmp1, \acc1, \tmp1

        or      \tmp1, \tmp1, \tmp2
        add     \acc2, \acc2, \tmp1
    .endm

    .globl fn_mul_4x4_u512_ct
    .type  fn_mul_4x4_u512_ct, @function

fn_mul_4x4_u512_ct:
    # save only s0..s3
    addi    sp, sp, -32
    sd      s0, 24(sp)
    sd      s1, 16(sp)
    sd      s2, 8(sp)
    sd      s3, 0(sp)

    # move out ptr away from a0 so a0 can be reused as temp
    mv      t1, a0

    # load a -> a4..a7
    ld      a4,  0(a1)     # a[0]
    ld      a5,  8(a1)     # a[1]
    ld      a6, 16(a1)     # a[2]
    ld      a7, 24(a1)     # a[3]

    # load b -> s0..s3
    ld      s0,  0(a2)     # b[0]
    ld      s1,  8(a2)     # b[1]
    ld      s2, 16(a2)     # b[2]
    ld      s3, 24(a2)     # b[3]

    # current column accumulator history
    li      t4, 0          # acc0
    li      t5, 0          # acc1

    # =============================================================
    # k = 0 : a0*b0
    # =============================================================
    li      t6, 0

    mul     t2, a4, s0
    mulhu   t3, a4, s0
    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a1

    sd      t4, 0(t1)
    mv      t4, t5
    mv      t5, t6

    # =============================================================
    # k = 1 : a0*b1 + a1*b0
    # =============================================================
    li      t6, 0

    mul     t2, a4, s1
    mulhu   t3, a4, s1
    mul     a0, a5, s0         # precompute next lo
    mulhu   a1, a5, s0         # precompute next hi

    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2
    ACC192_ADD128 t4, t5, t6, a0, a1, t0, a2

    sd      t4, 8(t1)
    mv      t4, t5
    mv      t5, t6

    # =============================================================
    # k = 2 : a0*b2 + a1*b1 + a2*b0
    # schedule: compute pair1 + pair2 first, acc pair1, compute pair3, acc pair2, acc pair3
    # =============================================================
    li      t6, 0

    mul     t2, a4, s2         # p0.lo
    mulhu   t3, a4, s2         # p0.hi
    mul     a0, a5, s1         # p1.lo
    mulhu   a1, a5, s1         # p1.hi

    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2

    mul     t2, a6, s0         # p2.lo
    mulhu   t3, a6, s0         # p2.hi

    ACC192_ADD128 t4, t5, t6, a0, a1, t0, a2
    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2

    sd      t4, 16(t1)
    mv      t4, t5
    mv      t5, t6

    # =============================================================
    # k = 3 : a0*b3 + a1*b2 + a2*b1 + a3*b0
    # schedule: issue p0,p1; acc p0; issue p2; acc p1; issue p3; acc p2; acc p3
    # =============================================================
    li      t6, 0

    mul     t2, a4, s3         # p0.lo
    mulhu   t3, a4, s3         # p0.hi
    mul     a0, a5, s2         # p1.lo
    mulhu   a1, a5, s2         # p1.hi

    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2

    mul     t2, a6, s1         # p2.lo
    mulhu   t3, a6, s1         # p2.hi

    ACC192_ADD128 t4, t5, t6, a0, a1, t0, a2

    mul     a0, a7, s0         # p3.lo
    mulhu   a1, a7, s0         # p3.hi

    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2
    ACC192_ADD128 t4, t5, t6, a0, a1, t0, a2

    sd      t4, 24(t1)
    mv      t4, t5
    mv      t5, t6

    # =============================================================
    # k = 4 : a1*b3 + a2*b2 + a3*b1
    # =============================================================
    li      t6, 0

    mul     t2, a5, s3         # p0.lo
    mulhu   t3, a5, s3         # p0.hi
    mul     a0, a6, s2         # p1.lo
    mulhu   a1, a6, s2         # p1.hi

    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2

    mul     t2, a7, s1         # p2.lo
    mulhu   t3, a7, s1         # p2.hi

    ACC192_ADD128 t4, t5, t6, a0, a1, t0, a2
    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2

    sd      t4, 32(t1)
    mv      t4, t5
    mv      t5, t6

    # =============================================================
    # k = 5 : a2*b3 + a3*b2
    # =============================================================
    li      t6, 0

    mul     t2, a6, s3
    mulhu   t3, a6, s3
    mul     a0, a7, s2
    mulhu   a1, a7, s2

    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2
    ACC192_ADD128 t4, t5, t6, a0, a1, t0, a2

    sd      t4, 40(t1)
    mv      t4, t5
    mv      t5, t6

    # =============================================================
    # k = 6 : a3*b3
    # =============================================================
    li      t6, 0

    mul     t2, a7, s3
    mulhu   t3, a7, s3
    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a1

    sd      t4, 48(t1)
    mv      t4, t5
    mv      t5, t6

    # =============================================================
    # k = 7
    # =============================================================
    sd      t4, 56(t1)

    # restore s0..s3
    ld      s0, 24(sp)
    ld      s1, 16(sp)
    ld      s2, 8(sp)
    ld      s3, 0(sp)
    addi    sp, sp, 32
    ret

    .size fn_mul_4x4_u512_ct, .-fn_mul_4x4_u512_ct
