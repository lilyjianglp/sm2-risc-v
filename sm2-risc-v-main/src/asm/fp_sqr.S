# asm/fp_sqr.S
.text
.align 2

# -----------------------------------------------------------------------------
# (acc2:acc1:acc0) += lo + (hi<<64)
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# (acc2:acc1:acc0) += 2*(hi:lo)
# -----------------------------------------------------------------------------
.macro ACC192_ADD128_DBL acc0, acc1, acc2, lo, hi, tA, tB, tC
    srli    \tC, \hi, 63
    slli    \tA, \lo, 1
    slli    \tB, \hi, 1
    srli    \lo, \lo, 63
    or      \tB, \tB, \lo
    ACC192_ADD128 \acc0, \acc1, \acc2, \tA, \tB, \tA, \tB
    add     \acc2, \acc2, \tC
.endm


.globl fp_sqr
.type fp_sqr, @function

# -----------------------------------------------------------------------------
# void fp_sqr(fp_t* r, const fp_t* a)
#
# 寄存器轮转：
#   k=0  acc=(t6, t0, t1)
#   k=1  acc=(t0, t1, t6)   t6已store可复用
#   k=2  acc=(t1, t6, t0)   t0已store可复用
#   k=3  acc=(t6, t0, t1)   t1已store可复用  (3轮一循环)
#   k=4  acc=(t0, t1, t6)
#   k=5  acc=(t1, t6, t0)
#   k=6  acc=(t6, t0, t1)
#   k=7  t0 即最终进位，直接store
# -----------------------------------------------------------------------------
fp_sqr:
    addi    sp, sp, -128
    sd      ra, 120(sp)
    sd      s0, 112(sp)
    sd      s1, 104(sp)
    sd      s2,  96(sp)
    sd      s3,  88(sp)
    sd      s11, 80(sp)

    mv      s11, a0

    ld      s0,  0(a1)
    ld      s1,  8(a1)
    ld      s2, 16(a1)
    ld      s3, 24(a1)

    # ---------------- k=0  acc=(t6, t0, t1) ----------------
    li      t6, 0
    li      t0, 0
    li      t1, 0

    mul     t2, s0, s0
    mulhu   t3, s0, s0
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 0(sp)

    # ---------------- k=1  acc=(t0, t1, t6) ----------------
    li      t6, 0

    mul     t2, s0, s1
    mulhu   t3, s0, s1
    ACC192_ADD128_DBL t0, t1, t6, t2, t3, a3, a4, a5

    sd      t0, 8(sp)

    # ---------------- k=2  acc=(t1, t6, t0) ----------------
    li      t0, 0

    mul     t2, s0, s2
    mulhu   t3, s0, s2
    ACC192_ADD128_DBL t1, t6, t0, t2, t3, a3, a4, a5

    mul     t2, s1, s1
    mulhu   t3, s1, s1
    ACC192_ADD128 t1, t6, t0, t2, t3, a3, a4

    sd      t1, 16(sp)

    # ---------------- k=3  acc=(t6, t0, t1) ----------------
    li      t1, 0

    mul     t2, s0, s3
    mulhu   t3, s0, s3
    ACC192_ADD128_DBL t6, t0, t1, t2, t3, a3, a4, a5

    mul     t2, s1, s2
    mulhu   t3, s1, s2
    ACC192_ADD128_DBL t6, t0, t1, t2, t3, a3, a4, a5

    sd      t6, 24(sp)

    # ---------------- k=4  acc=(t0, t1, t6) ----------------
    li      t6, 0

    mul     t2, s1, s3
    mulhu   t3, s1, s3
    ACC192_ADD128_DBL t0, t1, t6, t2, t3, a3, a4, a5

    mul     t2, s2, s2
    mulhu   t3, s2, s2
    ACC192_ADD128 t0, t1, t6, t2, t3, a3, a4

    sd      t0, 32(sp)

    # ---------------- k=5  acc=(t1, t6, t0) ----------------
    li      t0, 0

    mul     t2, s2, s3
    mulhu   t3, s2, s3
    ACC192_ADD128_DBL t1, t6, t0, t2, t3, a3, a4, a5

    sd      t1, 40(sp)

    # ---------------- k=6  acc=(t6, t0, t1) ----------------
    li      t1, 0

    mul     t2, s3, s3
    mulhu   t3, s3, s3
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 48(sp)

    # ---------------- k=7  只剩进位 ----------------
    sd      t0, 56(sp)

    # call fp_reduce(r, T)
    mv      a0, s11
    mv      a1, sp
    call    fp_reduce

    ld      ra, 120(sp)
    ld      s0, 112(sp)
    ld      s1, 104(sp)
    ld      s2,  96(sp)
    ld      s3,  88(sp)
    ld      s11, 80(sp)
    addi    sp, sp, 128
    ret
.size fp_sqr, .-fp_sqr