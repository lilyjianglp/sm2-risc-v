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
    sltu    \tmp1, \acc0, \lo      # tmp1 = carry0 (from low add)

    # acc1 += hi
    add     \acc1, \acc1, \hi
    sltu    \tmp2, \acc1, \hi      # tmp2 = carry1 (from high add)

    # acc1 += carry0
    add     \acc1, \acc1, \tmp1
    sltu    \tmp1, \acc1, \tmp1    # tmp1 = carry2 (from carry add)

    # acc2 += carry1 | carry2
    or      \tmp1, \tmp1, \tmp2
    add     \acc2, \acc2, \tmp1
.endm

.text
.align 2
.globl fp_mul
.type fp_mul, @function

# -----------------------------------------------------------------------------
# void fp_mul(fp_t* r, const fp_t* a, const fp_t* b)
#
# ABI:
#   a0 = r (结果地址)
#   a1 = a (操作数1地址)
#   a2 = b (操作数2地址)
# -----------------------------------------------------------------------------
fp_mul:
    # 1. 栈帧调整与寄存器保存
    # 需要保存: ra, s0-s11 (13个)
    # 需要空间: T[0]..T[7] (8个uint64 = 64 bytes)
    # 总共: 13*8 + 64 = 104 + 64 = 168 bytes
    # 对齐到 16 字节 -> 176 bytes
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
    sd      s8, 96(sp)
    sd      s9, 88(sp)
    sd      s10, 80(sp)
    sd      s11, 72(sp)
    
    # 栈布局:
    # 0(sp) - 63(sp): T[0]..T[7] (中间乘积)
    
    # 保存结果指针 r (a0) 到 s11，因为 a0 稍后会被覆盖，且最后调用 fp_reduce 需要用
    mv      s11, a0

    # 2. 加载输入
    # load a -> s0..s3
    ld      s0,  0(a1)
    ld      s1,  8(a1)
    ld      s2, 16(a1)
    ld      s3, 24(a1)

    # load b -> s4..s7
    ld      s4,  0(a2)
    ld      s5,  8(a2)
    ld      s6, 16(a2)
    ld      s7, 24(a2)

    # 3. 乘法阶段 (Comba Method) - 计算 512 位乘积
    
    # 初始状态
    li      t4, 0      # c0 (acc1 prev)
    li      t5, 0      # c1 (acc2 prev)

    # --- k = 0 ---
    mv      t6, t4     # acc0
    mv      t0, t5     # acc1
    li      t1, 0      # acc2

    # a0 * b0
    mul     t2, s0, s4
    mulhu   t3, s0, s4
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 0(sp)     # Store T[0]

    mv      t4, t0        # update c0
    mv      t5, t1        # update c1

    # --- k = 1 ---
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

    # a0 * b1
    mul     t2, s0, s5
    mulhu   t3, s0, s5
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    # a1 * b0
    mul     t2, s1, s4
    mulhu   t3, s1, s4
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 8(sp)     # Store T[1]

    mv      t4, t0
    mv      t5, t1

    # --- k = 2 ---
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

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

    sd      t6, 16(sp)    # Store T[2]

    mv      t4, t0
    mv      t5, t1

    # --- k = 3 ---
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

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

    sd      t6, 24(sp)    # Store T[3]

    mv      t4, t0
    mv      t5, t1

    # --- k = 4 ---
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

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

    sd      t6, 32(sp)    # Store T[4]

    mv      t4, t0
    mv      t5, t1

    # --- k = 5 ---
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

    # a2 * b3
    mul     t2, s2, s7
    mulhu   t3, s2, s7
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    # a3 * b2
    mul     t2, s3, s6
    mulhu   t3, s3, s6
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 40(sp)    # Store T[5]

    mv      t4, t0
    mv      t5, t1

    # --- k = 6 ---
    mv      t6, t4
    mv      t0, t5
    li      t1, 0

    # a3 * b3
    mul     t2, s3, s7
    mulhu   t3, s3, s7
    ACC192_ADD128 t6, t0, t1, t2, t3, a3, a4

    sd      t6, 48(sp)    # Store T[6]

    mv      t4, t0
    mv      t5, t1

    # --- k = 7 ---
    # 此时只剩下进位 c0
    sd      t4, 56(sp)    # Store T[7]

    # 4. 调用约减函数 fp_reduce
    # ABI: fp_reduce(fp_t* r, uint64_t* T)
    # a0: 结果地址 r (之前保存在 s11)
    # a1: 中间结果 T 的地址 (sp)
    
    mv      a0, s11       # 恢复 r 的地址
    mv      a1, sp        # T 的地址就是当前栈指针指向的位置
    call    fp_reduce

    # 5. 恢复现场并返回
    ld      ra, 168(sp)
    ld      s0, 160(sp)
    ld      s1, 152(sp)
    ld      s2, 144(sp)
    ld      s3, 136(sp)
    ld      s4, 128(sp)
    ld      s5, 120(sp)
    ld      s6, 112(sp)
    ld      s7, 104(sp)
    ld      s8, 96(sp)
    ld      s9, 88(sp)
    ld      s10, 80(sp)
    ld      s11, 72(sp)
    addi    sp, sp, 176
    ret
    .size fp_mul, .-fp_mul

   