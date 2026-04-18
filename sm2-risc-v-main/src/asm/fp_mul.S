# -----------------------------------------------------------------------------
# 宏定义：ACC192_ADD128
# 作用：(acc2:acc1:acc0) += lo + (hi << 64)
# 参数：
#   acc0, acc1, acc2: 累加器寄存器 (acc0=低, acc1=中, acc2=高)
#   lo, hi:           要加的 128 位数 (lo=低64, hi=高64)
#   tmp1, tmp2:       临时寄存器 (破坏性使用)
# -----------------------------------------------------------------------------
.macro ACC192_ADD128 acc0, acc1, acc2, lo, hi, tmp1, tmp2
    # acc0 += lo
    add     \acc0, \acc0, \lo
    sltu    \tmp1, \acc0, \lo      # tmp1 = carry0

    # acc1 += hi
    add     \acc1, \acc1, \hi
    sltu    \tmp2, \acc1, \hi      # tmp2 = carry1

    # acc1 += carry0
    add     \acc1, \acc1, \tmp1
    sltu    \tmp1, \acc1, \tmp1    # tmp1 = carry2

    # acc2 += carry1 | carry2
    or      \tmp1, \tmp1, \tmp2
    add     \acc2, \acc2, \tmp1
.endm

# -----------------------------------------------------------------------------
# 宏定义：ACC192_SHIFT
# 作用：Comba 列推进
#       当前列输出 acc0 后，下一列初始值变为：
#       acc0 <- acc1
#       acc1 <- acc2
#       acc2 <- 0
# -----------------------------------------------------------------------------
.macro ACC192_SHIFT acc0, acc1, acc2
    mv      \acc0, \acc1
    mv      \acc1, \acc2
    li      \acc2, 0
.endm

.text
.align  2
.globl  fp_mul
.type   fp_mul, @function

# -----------------------------------------------------------------------------
# void fp_mul(fp_t* r, const fp_t* a, const fp_t* b)
#
# ABI:
#   a0 = r
#   a1 = a
#   a2 = b
# -----------------------------------------------------------------------------
fp_mul:
    # 1. 建栈并保存现场
    # 保存: ra, s0-s11 = 13 * 8 = 104 bytes
    # T[0..7] = 8 * 8 = 64 bytes
    # 总计 168 bytes，向上按 16 对齐到 176 bytes
    addi    sp, sp, -176
    sd      ra, 168(sp)
    sd      s0, 160(sp)
    sd      s1, 152(sp)
    sd      s2, 144(sp)
    sd      s3, 136(sp)
    sd      s4, 128(sp)
    sd      s5, 120(sp)
    sd      s6, 112(sp)
    sd      s7, 104(sp)
    sd      s8,  96(sp)
    sd      s9,  88(sp)
    sd      s10, 80(sp)
    sd      s11, 72(sp)

    # 栈布局:
    #   0(sp)  ~ 56(sp): T[0]..T[7]

    # 保存结果指针 r
    mv      s11, a0

    # 2. 加载 a, b
    ld      s0,  0(a1)
    ld      s1,  8(a1)
    ld      s2, 16(a1)
    ld      s3, 24(a1)

    ld      s4,  0(a2)
    ld      s5,  8(a2)
    ld      s6, 16(a2)
    ld      s7, 24(a2)

    # 3. Comba 乘法
    # 固定:
    #   t6 = acc0
    #   t0 = acc1
    #   t1 = acc2
    #   t2 = lo
    #   t3 = hi
    #   a3, a4 = carry 临时寄存器
    li      t6, 0
    li      t0, 0
    li      t1, 0

    # --- k = 0 ---
    # a0 * b0
    mul     t2, s0, s4
    mulhu   t3, s0, s4
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 0(sp)             # T[0]
    ACC192_SHIFT t6, t0, t1

    # --- k = 1 ---
    # a0 * b1
    mul     t2, s0, s5
    mulhu   t3, s0, s5
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    # a1 * b0
    mul     t2, s1, s4
    mulhu   t3, s1, s4
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 8(sp)             # T[1]
    ACC192_SHIFT t6, t0, t1

    # --- k = 2 ---
    # a0 * b2
    mul     t2, s0, s6
    mulhu   t3, s0, s6
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    # a1 * b1
    mul     t2, s1, s5
    mulhu   t3, s1, s5
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    # a2 * b0
    mul     t2, s2, s4
    mulhu   t3, s2, s4
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 16(sp)            # T[2]
    ACC192_SHIFT t6, t0, t1

    # --- k = 3 ---
    # a0 * b3
    mul     t2, s0, s7
    mulhu   t3, s0, s7
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    # a1 * b2
    mul     t2, s1, s6
    mulhu   t3, s1, s6
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    # a2 * b1
    mul     t2, s2, s5
    mulhu   t3, s2, s5
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    # a3 * b0
    mul     t2, s3, s4
    mulhu   t3, s3, s4
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 24(sp)            # T[3]
    ACC192_SHIFT t6, t0, t1

    # --- k = 4 ---
    # a1 * b3
    mul     t2, s1, s7
    mulhu   t3, s1, s7
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    # a2 * b2
    mul     t2, s2, s6
    mulhu   t3, s2, s6
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    # a3 * b1
    mul     t2, s3, s5
    mulhu   t3, s3, s5
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 32(sp)            # T[4]
    ACC192_SHIFT t6, t0, t1

    # --- k = 5 ---
    # a2 * b3
    mul     t2, s2, s7
    mulhu   t3, s2, s7
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    # a3 * b2
    mul     t2, s3, s6
    mulhu   t3, s3, s6
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 40(sp)            # T[5]
    ACC192_SHIFT t6, t0, t1

    # --- k = 6 ---
    # a3 * b3
    mul     t2, s3, s7
    mulhu   t3, s3, s7
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 48(sp)            # T[6]
    ACC192_SHIFT t6, t0, t1

    # --- k = 7 ---
    # 最后一列只剩 1 个 limb
    sd      t6, 56(sp)            # T[7]

    # 4. 调用约减函数 fp_reduce(r, T)
    mv      a0, s11
    mv      a1, sp
    call    fp_reduce

    # 5. 恢复现场
    ld      ra, 168(sp)
    ld      s0, 160(sp)
    ld      s1, 152(sp)
    ld      s2, 144(sp)
    ld      s3, 136(sp)
    ld      s4, 128(sp)
    ld      s5, 120(sp)
    ld      s6, 112(sp)
    ld      s7, 104(sp)
    ld      s8,  96(sp)
    ld      s9,  88(sp)
    ld      s10, 80(sp)
    ld      s11, 72(sp)
    addi    sp, sp, 176
    ret

.size fp_mul, .-fp_mul