# =============================================================================
# fn_mont.S
#
# RISC-V RV64 + M extension
#
# Implements:
#   void fn_mont_mul_ct(fn_t *r, const fn_t *a, const fn_t *b);
#
# Notes:
# - merged kernel:
#     1) 4x4 limb multiply -> t[0..7]
#     2) t[8] = 0
#     3) Montgomery reduction for SM2 order n
#
# ABI:
#   a0 = r
#   a1 = a
#   a2 = b
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


# =============================================================================
# void fn_mont_mul_ct(fn_t *r, const fn_t *a, const fn_t *b)
#
# Stack layout (128 bytes total):
#   sp+0   .. sp+71    : t[9] (9 * 8 = 72 bytes)
#   sp+80              : saved r pointer
#   sp+88  .. sp+120   : saved s0..s4
# =============================================================================
    .globl fn_mont_mul_ct
    .type  fn_mont_mul_ct, @function
fn_mont_mul_ct:
    # reserve stack
    addi    sp, sp, -128

    # save r pointer + callee-saved regs
    sd      a0,  80(sp)
    sd      s0,  88(sp)
    sd      s1,  96(sp)
    sd      s2, 104(sp)
    sd      s3, 112(sp)
    sd      s4, 120(sp)

    # -------------------------------------------------------------------------
    # Part 1: 4x4 -> 8-limb multiply
    #   t[0..7] = a * b
    #   t[8]    = 0
    # -------------------------------------------------------------------------

    # t base on stack
    mv      t1, sp                  # t1 = &t[0]

    # load a -> a4..a7
    ld      a4,  0(a1)              # a[0]
    ld      a5,  8(a1)              # a[1]
    ld      a6, 16(a1)              # a[2]
    ld      a7, 24(a1)              # a[3]

    # load b -> s0..s3
    ld      s0,  0(a2)              # b[0]
    ld      s1,  8(a2)              # b[1]
    ld      s2, 16(a2)              # b[2]
    ld      s3, 24(a2)              # b[3]

    # rolling accumulator history
    li      t4, 0                   # acc0
    li      t5, 0                   # acc1

    # -------------------------------------------------------------------------
    # k = 0 : a0*b0
    # -------------------------------------------------------------------------
    li      t6, 0                   # acc2

    mul     t2, a4, s0
    mulhu   t3, a4, s0
    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a1

    sd      t4, 0(t1)
    mv      t4, t5
    mv      t5, t6

    # -------------------------------------------------------------------------
    # k = 1 : a0*b1 + a1*b0
    # -------------------------------------------------------------------------
    li      t6, 0

    mul     t2, a4, s1
    mulhu   t3, a4, s1
    mul     a0, a5, s0
    mulhu   a1, a5, s0

    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2
    ACC192_ADD128 t4, t5, t6, a0, a1, t0, a2

    sd      t4, 8(t1)
    mv      t4, t5
    mv      t5, t6

    # -------------------------------------------------------------------------
    # k = 2 : a0*b2 + a1*b1 + a2*b0
    # -------------------------------------------------------------------------
    li      t6, 0

    mul     t2, a4, s2
    mulhu   t3, a4, s2
    mul     a0, a5, s1
    mulhu   a1, a5, s1

    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2

    mul     t2, a6, s0
    mulhu   t3, a6, s0

    ACC192_ADD128 t4, t5, t6, a0, a1, t0, a2
    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2

    sd      t4, 16(t1)
    mv      t4, t5
    mv      t5, t6

    # -------------------------------------------------------------------------
    # k = 3 : a0*b3 + a1*b2 + a2*b1 + a3*b0
    # -------------------------------------------------------------------------
    li      t6, 0

    mul     t2, a4, s3
    mulhu   t3, a4, s3
    mul     a0, a5, s2
    mulhu   a1, a5, s2

    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2

    mul     t2, a6, s1
    mulhu   t3, a6, s1

    ACC192_ADD128 t4, t5, t6, a0, a1, t0, a2

    mul     a0, a7, s0
    mulhu   a1, a7, s0

    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2
    ACC192_ADD128 t4, t5, t6, a0, a1, t0, a2

    sd      t4, 24(t1)
    mv      t4, t5
    mv      t5, t6

    # -------------------------------------------------------------------------
    # k = 4 : a1*b3 + a2*b2 + a3*b1
    # -------------------------------------------------------------------------
    li      t6, 0

    mul     t2, a5, s3
    mulhu   t3, a5, s3
    mul     a0, a6, s2
    mulhu   a1, a6, s2

    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2

    mul     t2, a7, s1
    mulhu   t3, a7, s1

    ACC192_ADD128 t4, t5, t6, a0, a1, t0, a2
    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a2

    sd      t4, 32(t1)
    mv      t4, t5
    mv      t5, t6

    # -------------------------------------------------------------------------
    # k = 5 : a2*b3 + a3*b2
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # k = 6 : a3*b3
    # -------------------------------------------------------------------------
    li      t6, 0

    mul     t2, a7, s3
    mulhu   t3, a7, s3
    ACC192_ADD128 t4, t5, t6, t2, t3, t0, a1

    sd      t4, 48(t1)
    mv      t4, t5
    mv      t5, t6

    # -------------------------------------------------------------------------
    # k = 7
    # -------------------------------------------------------------------------
    sd      t4, 56(t1)
    sd      zero, 64(t1)           # t[8] = 0

    # -------------------------------------------------------------------------
    # Part 2: Montgomery reduction
    #   restore:
    #       a0 = r
    #       a1 = &t[0]
    # -------------------------------------------------------------------------
    ld      a0, 80(sp)             # restore r pointer
    mv      a1, sp                 # a1 = &t[0]

    # load constants
    la      t6, .L_n
    ld      s0,  0(t6)             # n0
    ld      s1,  8(t6)             # n1
    ld      s2, 16(t6)             # n2
    ld      s3, 24(t6)             # n3

    la      t0, .L_n0prime
    ld      s4, 0(t0)              # n0prime

    la      t5, .L_mont_one        # &mont_one[0]

    # for (i=0; i<4; i++)
    li      t0, 0
.L_i_loop_mont_mul:
    slli    t1, t0, 3
    add     t2, a1, t1             # t2 = &t[i]

    # m = (t[i] * n0prime) low64
    ld      a2, 0(t2)
    mul     a3, a2, s4             # m
    li      a4, 0                  # carry = 0

    # j = 0
    mul     a6, a3, s0
    mulhu   a4, a3, s0
    ld      a2, 0(t2)
    add     a6, a6, a2
    sltu    a2, a6, a2
    add     a4, a4, a2
    sd      a6, 0(t2)

    # j = 1
    mul     a6, a3, s1
    mulhu   a7, a3, s1
    ld      a2, 8(t2)
    add     a6, a6, a2
    sltu    a2, a6, a2
    add     a7, a7, a2
    add     a6, a6, a4
    sltu    a2, a6, a4
    add     a4, a7, a2
    sd      a6, 8(t2)

    # j = 2
    mul     a6, a3, s2
    mulhu   a7, a3, s2
    ld      a2, 16(t2)
    add     a6, a6, a2
    sltu    a2, a6, a2
    add     a7, a7, a2
    add     a6, a6, a4
    sltu    a2, a6, a4
    add     a4, a7, a2
    sd      a6, 16(t2)

    # j = 3
    mul     a6, a3, s3
    mulhu   a7, a3, s3
    ld      a2, 24(t2)
    add     a6, a6, a2
    sltu    a2, a6, a2
    add     a7, a7, a2
    add     a6, a6, a4
    sltu    a2, a6, a4
    add     a4, a7, a2
    sd      a6, 24(t2)

    # fixed propagation to t[i+4..8]
    ld      a2, 32(t2)
    add     a3, a2, a4
    sltu    a4, a3, a2
    sd      a3, 32(t2)

    ld      a2, 40(t2)
    add     a3, a2, a4
    sltu    a4, a3, a2
    sd      a3, 40(t2)

    ld      a2, 48(t2)
    add     a3, a2, a4
    sltu    a4, a3, a2
    sd      a3, 48(t2)

    ld      a2, 56(t2)
    add     a3, a2, a4
    sltu    a4, a3, a2
    sd      a3, 56(t2)

    ld      a2, 64(t2)
    add     a3, a2, a4
    sltu    a4, a3, a2
    sd      a3, 64(t2)

    addi    t0, t0, 1
    li      t1, 4
    blt     t0, t1, .L_i_loop_mont_mul

    # r = t[4..7]
    ld      t0, 32(a1)
    sd      t0, 0(a0)
    ld      t0, 40(a1)
    sd      t0, 8(a0)
    ld      t0, 48(a1)
    sd      t0, 16(a0)
    ld      t0, 56(a1)
    sd      t0, 24(a0)

    # t8 mask
    ld      t0, 64(a1)
    andi    t0, t0, 1
    neg     t0, t0

    # load r
    ld      t1, 0(a0)
    ld      t2, 8(a0)
    ld      t3, 16(a0)
    ld      t4, 24(a0)

    # add mont_one & mask
    ld      a2, 0(t5)
    and     a2, a2, t0
    add     t1, t1, a2
    sltu    a7, t1, a2

    ld      a2, 8(t5)
    and     a2, a2, t0
    add     t2, t2, a2
    sltu    a6, t2, a2
    add     t2, t2, a7
    sltu    a7, t2, a7
    or      a7, a7, a6

    # mont_one[2] == 0, so only propagate carry
    add     t3, t3, a7
    sltu    a7, t3, a7

    ld      a2, 24(t5)
    and     a2, a2, t0
    add     t4, t4, a2
    sltu    a6, t4, a2
    add     t4, t4, a7
    sltu    a7, t4, a7
    or      a7, a7, a6

    sd      t1, 0(a0)
    sd      t2, 8(a0)
    sd      t3, 16(a0)
    sd      t4, 24(a0)

    # -------- cond_sub_n #1 --------
    ld      t1, 0(a0)
    ld      t2, 8(a0)
    ld      t3, 16(a0)
    ld      t4, 24(a0)

    mv      a2, t1
    mv      a3, t2
    mv      a4, t3
    mv      a5, t4

    sub     t1, t1, s0
    sltu    t0, a2, s0

    sub     t2, t2, s1
    sltu    a6, a3, s1
    sub     t2, t2, t0
    sltu    t0, t2, t0
    or      t0, t0, a6

    sub     t3, t3, s2
    sltu    a6, a4, s2
    sub     t3, t3, t0
    sltu    t0, t3, t0
    or      t0, t0, a6

    sub     t4, t4, s3
    sltu    a6, a5, s3
    sub     t4, t4, t0
    sltu    t0, t4, t0
    or      t0, t0, a6

    xori    t0, t0, 1
    neg     t0, t0
    not     a6, t0

    and     t1, t1, t0
    and     a2, a2, a6
    or      t1, t1, a2
    sd      t1, 0(a0)

    and     t2, t2, t0
    and     a3, a3, a6
    or      t2, t2, a3
    sd      t2, 8(a0)

    and     t3, t3, t0
    and     a4, a4, a6
    or      t3, t3, a4
    sd      t3, 16(a0)

    and     t4, t4, t0
    and     a5, a5, a6
    or      t4, t4, a5
    sd      t4, 24(a0)

    # -------- cond_sub_n #2 --------
    ld      t1, 0(a0)
    ld      t2, 8(a0)
    ld      t3, 16(a0)
    ld      t4, 24(a0)

    mv      a2, t1
    mv      a3, t2
    mv      a4, t3
    mv      a5, t4

    sub     t1, t1, s0
    sltu    t0, a2, s0

    sub     t2, t2, s1
    sltu    a6, a3, s1
    sub     t2, t2, t0
    sltu    t0, t2, t0
    or      t0, t0, a6

    sub     t3, t3, s2
    sltu    a6, a4, s2
    sub     t3, t3, t0
    sltu    t0, t3, t0
    or      t0, t0, a6

    sub     t4, t4, s3
    sltu    a6, a5, s3
    sub     t4, t4, t0
    sltu    t0, t4, t0
    or      t0, t0, a6

    xori    t0, t0, 1
    neg     t0, t0
    not     a6, t0

    and     t1, t1, t0
    and     a2, a2, a6
    or      t1, t1, a2
    sd      t1, 0(a0)

    and     t2, t2, t0
    and     a3, a3, a6
    or      t2, t2, a3
    sd      t2, 8(a0)

    and     t3, t3, t0
    and     a4, a4, a6
    or      t3, t3, a4
    sd      t3, 16(a0)

    and     t4, t4, t0
    and     a5, a5, a6
    or      t4, t4, a5
    sd      t4, 24(a0)

    # restore
    ld      s0,  88(sp)
    ld      s1,  96(sp)
    ld      s2, 104(sp)
    ld      s3, 112(sp)
    ld      s4, 120(sp)
    addi    sp, sp, 128
    ret

    .size fn_mont_mul_ct, .-fn_mont_mul_ct


# =============================================================================
# constants
# =============================================================================
    .section .rodata
    .align 3
.L_n:
    .quad   0x53BBF40939D54123
    .quad   0x7203DF6B21C6052B
    .quad   0xFFFFFFFFFFFFFFFF
    .quad   0xFFFFFFFEFFFFFFFF

.L_n0prime:
    .quad   0x327f9e8872350975

.L_mont_one:
    .quad   0xAC440BF6C62ABEDD
    .quad   0x8DFC2094DE39FAD4
    .quad   0x0000000000000000
    .quad   0x0000000100000000
