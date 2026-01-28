# asm/fp_sqr.S
.text
.align 2

# -----------------------------------------------------------------------------
# (acc2:acc1:acc0) += lo + (hi<<64)
# -----------------------------------------------------------------------------
.macro ACC192_ADD128 acc0, acc1, acc2, lo, hi, tmp1, tmp2
    add     \acc0, \acc0, \lo
    sltu    \tmp1, \acc0, \lo      # carry0

    add     \acc1, \acc1, \hi
    sltu    \tmp2, \acc1, \hi      # carry1

    add     \acc1, \acc1, \tmp1
    sltu    \tmp1, \acc1, \tmp1    # carry2 (from adding carry0)

    or      \tmp1, \tmp1, \tmp2
    add     \acc2, \acc2, \tmp1
.endm

# -----------------------------------------------------------------------------
# (acc2:acc1:acc0) += 2*(hi:lo)
# tmp regs: tA,tB,tC used
# -----------------------------------------------------------------------------
.macro ACC192_ADD128_DBL acc0, acc1, acc2, lo, hi, tA, tB, tC
    # carry128 = hi >> 63
    srli    \tC, \hi, 63

    # dlo = lo << 1
    slli    \tA, \lo, 1

    # dhi = (hi<<1) | (lo>>63)
    slli    \tB, \hi, 1
    srli    \lo, \lo, 63           # reuse lo as temp: lo>>63
    or      \tB, \tB, \lo

    # add (dhi:dlo)
    ACC192_ADD128 \acc0, \acc1, \acc2, \tA, \tB, \tA, \tB

    # acc2 += carry128
    add     \acc2, \acc2, \tC
.endm


.globl fp_sqr
.type fp_sqr, @function
# -----------------------------------------------------------------------------
# void fp_sqr(fp_t* r, const fp_t* a)
# ABI:
#   a0 = r
#   a1 = a
# -----------------------------------------------------------------------------
fp_sqr:
    # stack:
    #   0..63   : T[0..7]
    #   64..    : saved regs
    # save ra, s0-s3, s11  => 6*8 = 48
    # total 64 + 48 = 112, align to 16 => 128
    addi    sp, sp, -128
    sd      ra, 120(sp)
    sd      s0, 112(sp)
    sd      s1, 104(sp)
    sd      s2, 96(sp)
    sd      s3, 88(sp)
    sd      s11, 80(sp)

    mv      s11, a0          # keep r

    # load a0..a3 into s0..s3
    ld      s0,  0(a1)
    ld      s1,  8(a1)
    ld      s2, 16(a1)
    ld      s3, 24(a1)

    # c0=t4, c1=t5
    li      t4, 0
    li      t5, 0

    # ---------------- k = 0 : a0*a0 ----------------
    mv      t6, t4           # acc0
    mv      t0, t5           # acc1
    li      t1, 0            # acc2

    mul     t2, s0, s0
    mulhu   t3, s0, s0
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 0(sp)
    mv      t4, t0
    mv      t5, t1

    # ---------------- k = 1 : 2*a0*a1 ----------------
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

    mul     t2, s0, s1
    mulhu   t3, s0, s1
    # dbl add
    # NOTE: ACC192_ADD128_DBL will clobber t2/t3 temps internally; we pass lo=t2 hi=t3
    ACC192_ADD128_DBL t6, t0, t1, t2, t3, a3, a4, a5

    sd      t6, 8(sp)
    mv      t4, t0
    mv      t5, t1

    # ---------------- k = 2 : 2*a0*a2 + a1*a1 ----------------
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

    mul     t2, s0, s2
    mulhu   t3, s0, s2
    ACC192_ADD128_DBL t6, t0, t1, t2, t3, a3, a4, a5

    mul     t2, s1, s1
    mulhu   t3, s1, s1
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 16(sp)
    mv      t4, t0
    mv      t5, t1

    # ---------------- k = 3 : 2*a0*a3 + 2*a1*a2 ----------------
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

    mul     t2, s0, s3
    mulhu   t3, s0, s3
    ACC192_ADD128_DBL t6, t0, t1, t2, t3, a3, a4, a5

    mul     t2, s1, s2
    mulhu   t3, s1, s2
    ACC192_ADD128_DBL t6, t0, t1, t2, t3, a3, a4, a5

    sd      t6, 24(sp)
    mv      t4, t0
    mv      t5, t1

    # ---------------- k = 4 : 2*a1*a3 + a2*a2 ----------------
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

    mul     t2, s1, s3
    mulhu   t3, s1, s3
    ACC192_ADD128_DBL t6, t0, t1, t2, t3, a3, a4, a5

    mul     t2, s2, s2
    mulhu   t3, s2, s2
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 32(sp)
    mv      t4, t0
    mv      t5, t1

    # ---------------- k = 5 : 2*a2*a3 ----------------
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

    mul     t2, s2, s3
    mulhu   t3, s2, s3
    ACC192_ADD128_DBL t6, t0, t1, t2, t3, a3, a4, a5

    sd      t6, 40(sp)
    mv      t4, t0
    mv      t5, t1

    # ---------------- k = 6 : a3*a3 ----------------
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

    mul     t2, s3, s3
    mulhu   t3, s3, s3
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 48(sp)
    mv      t4, t0
    mv      t5, t1

    # ---------------- k = 7 : carry ----------------
    sd      t4, 56(sp)

    # call fp_reduce(r, T)
    mv      a0, s11
    mv      a1, sp
    call    fp_reduce

    # restore
    ld      ra, 120(sp)
    ld      s0, 112(sp)
    ld      s1, 104(sp)
    ld      s2, 96(sp)
    ld      s3, 88(sp)
    ld      s11, 80(sp)
    addi    sp, sp, 128
    ret
.size fp_sqr, .-fp_sqr
