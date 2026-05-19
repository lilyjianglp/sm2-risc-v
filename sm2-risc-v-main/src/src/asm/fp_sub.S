     /* fp_sub */
    .globl  fp_sub
    .type   fp_sub, @function

# void fp_sub(fp_t* r=a0, const fp_t* a=a1, const fp_t* b=a2)
fp_sub:
    # ---------------------------------------------------------
    # Load a[0..3] into t0..t3
    # Load b[0..3] into t4..t6 and a7 (RISC-V has t0..t6 only)
    # ---------------------------------------------------------
    ld      t0,   0(a1)          # a0
    ld      t1,   8(a1)          # a1
    ld      t2,  16(a1)          # a2
    ld      t3,  24(a1)          # a3

    ld      t4,   0(a2)          # b0
    ld      t5,   8(a2)          # b1
    ld      t6,  16(a2)          # b2
    ld      a7,  24(a2)          # b3

    # ---------------------------------------------------------
    # 1) t = a - b  (borrow chain)
    #    t limbs -> a3,a4,a5,a6
    #    borrow  -> a2 (0/1)   (reuse a2 safely)
    #
    # subb_u64:
    #   diff  = a - b
    #   diff2 = diff - bin
    #   borrow = (a<b) | (diff<bin)   <-- IMPORTANT: diff before subtract bin
    # ---------------------------------------------------------
    sub     a3, t0, t4
    sltu    a2, t0, t4           # borrow0

    # limb1
    sub     t0, t1, t5           # diff  (keep in t0)
    sltu    t4, t1, t5           # borrow1 in t4
    sub     a4, t0, a2           # diff2 -> a4
    sltu    t5, t0, a2           # borrow2: (diff < bin)
    or      a2, t4, t5           # new borrow

    # limb2
    sub     t0, t2, t6
    sltu    t4, t2, t6
    sub     a5, t0, a2
    sltu    t5, t0, a2
    or      a2, t4, t5

    # limb3
    sub     t0, t3, a7
    sltu    t4, t3, a7
    sub     a6, t0, a2
    sltu    t5, t0, a2
    or      a2, t4, t5           # final borrow in a2

    # ---------------------------------------------------------
    # 2) Load p[0..3] from global FP_P into t4,t5,t6,a7
    # ---------------------------------------------------------
    la      t3, FP_P
    ld      t4,   0(t3)          # p0
    ld      t5,   8(t3)          # p1
    ld      t6,  16(t3)          # p2
    ld      a7,  24(t3)          # p3

    # ---------------------------------------------------------
    # 3) u = t + p  (carry chain)
    #    u limbs -> t0,t1,t2,t3
    #    carry   -> a1 (0/1)   (reuse a1 safely)
    # ---------------------------------------------------------
    add     t0, a3, t4
    sltu    a1, t0, a3           # carry0

    add     t1, a4, t5
    sltu    t4, t1, a4           # c_from_add
    add     t1, t1, a1
    sltu    t5, t1, a1           # c_from_cin
    or      a1, t4, t5           # carry1

    add     t2, a5, t6
    sltu    t4, t2, a5
    add     t2, t2, a1
    sltu    t5, t2, a1
    or      a1, t4, t5           # carry2

    add     t3, a6, a7
    add     t3, t3, a1           # final carry ignored (as in C)

    # ---------------------------------------------------------
    # 4) mask m = 0 - borrow
    #    m = ~0 if borrow==1 else 0
    #    ~m in a7
    # ---------------------------------------------------------
    neg     a2, a2               # a2 = m
    not     a7, a2               # a7 = ~m

    # ---------------------------------------------------------
    # 5) r = (u & m) | (t & ~m)   (branchless)
    # ---------------------------------------------------------
    and     t4, t0, a2
    and     t5, a3, a7
    or      t4, t4, t5
    sd      t4,   0(a0)

    and     t4, t1, a2
    and     t5, a4, a7
    or      t4, t4, t5
    sd      t4,   8(a0)

    and     t4, t2, a2
    and     t5, a5, a7
    or      t4, t4, t5
    sd      t4,  16(a0)

    and     t4, t3, a2
    and     t5, a6, a7
    or      t4, t4, t5
    sd      t4,  24(a0)

    ret

    .size fp_sub, .-fp_sub