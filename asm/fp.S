    .text
    .align  2
   
    /* fp_add */
    .globl  fp_add
    .type   fp_add, @function

# void fp_add(fp_t* r, const fp_t* a, const fp_t* b)
# a0 = r, a1 = a, a2 = b

fp_add:
    # ----------------------------
    # 1) Load a[0..3], b[0..3]
    # ----------------------------
    ld      a3,  0(a1)      # a0
    ld      a4,  8(a1)      # a1
    ld      a5, 16(a1)      # a2
    ld      a6, 24(a1)      # a3

    ld      t0,  0(a2)      # b0
    ld      t1,  8(a2)      # b1
    ld      t2, 16(a2)      # b2
    ld      t3, 24(a2)      # b3

    # a1/a2 pointers no longer needed -> can reuse a1/a2 as temporaries

    # ----------------------------
    # 2) t = a + b (with carry)
    # carry in t4 (0/1)
    # results: t0=t[0] in t5, t[1] in t6, t[2] in a2, t[3] in a7
    # ----------------------------
    li      t4, 0           # carry = 0

    # limb0: t5 = a3 + t0 + carry
    add     t5, a3, t0
    sltu    t4, t5, a3      # carry1
    

    # limb1: t6 = a4 + t1 + carry
    add     t6, a4, t1
    sltu    a2, t6, a4      # carry1
    mv      a1, t6
    add     t6, t6, t4
    sltu    a7, t6, a1      # carry2
    or      t4, a2, a7

    # limb2: a2 = a5 + t2 + carry
    add     a2, a5, t2
    sltu    a7, a2, a5      # carry1
    mv      a1, a2
    add     a2, a2, t4
    sltu    t1, a2, a1      # carry2  (t1 is free now)
    or      t4, a7, t1

    # limb3: a7 = a6 + t3 + carry
    add     a7, a6, t3
    sltu    t1, a7, a6      # carry1
    mv      a1, a7
    add     a7, a7, t4
    sltu    t2, a7, a1      # carry2  (t2 is free now)
    or      t4, t1, t2      # final carry in t4 (0/1)

    # Now:
    # t = (t5, t6, a2, a7)
    # carry = t4

    # ----------------------------
    # 3) u = t - p (with borrow)
    # borrow in t0 (0/1)
    # store u into a3,a4,a5,a6
    # p limbs:
    #   p0 = 0xFFFFFFFFFFFFFFFF
    #   p1 = 0xFFFFFFFF00000000
    #   p2 = 0xFFFFFFFFFFFFFFFF
    #   p3 = 0xFFFFFFFEFFFFFFFF
    # ----------------------------
    li      t0, 0           # borrow = 0

    # limb0: u0 = t5 - p0 - borrow
    li      t3, -1          # p0
    sub     a3, t5, t3      # diff = t5 - p0
    sltu    t0, t5, t3      # borrow1
    

    # limb1: u1 = t6 - p1 - borrow
    li      t3, 0xFFFFFFFF00000000
    sub     a1, t6, t3
    sltu    t2, t6, t3
    sub     a4, a1, t0
    sltu    t1, a1, t0
    or      t0, t2, t1

    # limb2: u2 = a2 - p2 - borrow
    li      t3, -1          # p2
    sub     a1, a2, t3
    sltu    t2, a2, t3
    sub     a5, a1, t0
    sltu    t1, a1, t0
    or      t0, t2, t1

    # limb3: u3 = a7 - p3 - borrow
    li      t3, 0xFFFFFFFEFFFFFFFF
    sub     a1, a7, t3
    sltu    t2, a7, t3
    sub     a6, a1, t0
    sltu    t1, a1, t0
    or      t0, t2, t1      # final borrow (0/1)

    # Now:
    # u = (a3,a4,a5,a6)
    # borrow = t0
    # carry  = t4

    # ----------------------------
    # 4) do_sub = carry | (1 ^ borrow)
    # mask m = 0 - do_sub  => 0xFFFF.. or 0x0
    # ----------------------------
    xori    t2, t0, 1       # (borrow ^ 1)
    or      t1, t4, t2      # do_sub (0/1)
    neg     t1, t1          # m = 0 - do_sub
    not     t2, t1          # ~m

    # ----------------------------
    # 5) r = (u & m) | (t & ~m)  (branchless)
    # store to *r
    # ----------------------------
    and     t3, a3, t1
    and     t4, t5, t2
    or      t3, t3, t4
    sd      t3,  0(a0)

    and     t3, a4, t1
    and     t4, t6, t2
    or      t3, t3, t4
    sd      t3,  8(a0)

    and     t3, a5, t1
    and     t4, a2, t2
    or      t3, t3, t4
    sd      t3, 16(a0)

    and     t3, a6, t1
    and     t4, a7, t2
    or      t3, t3, t4
    sd      t3, 24(a0)

    ret

    .size fp_add, .-fp_add
    
    
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

    
     
     
     
     /* fp_neg */
    .globl  fp_neg
    .type   fp_neg, @function

# void fp_neg(fp_t* r=a0, const fp_t* a=a1)
# r = (a==0) ? 0 : (p - a)
fp_neg:
    # ----------------------------
    # Load a[0..3] -> t0..t3
    # ----------------------------
    ld      t0,   0(a1)
    ld      t1,   8(a1)
    ld      t2,  16(a1)
    ld      t3,  24(a1)

    # ----------------------------
    # nz = (a != 0)
    # mask = 0 - nz   (all1 if nz else 0)
    # ----------------------------
    or      a2, t0, t1
    or      a2, a2, t2
    or      a2, a2, t3
    snez    a2, a2
    neg     a2, a2               # mask_nz in a2

    # ----------------------------
    # Load p[0..3] -> t4,t5,t6,a7
    # ----------------------------
    la      a6, FP_P
    ld      t4,   0(a6)          # p0
    ld      t5,   8(a6)          # p1
    ld      t6,  16(a6)          # p2
    ld      a7,  24(a6)          # p3

    # ----------------------------
    # t = p - a  (borrow chain, fp_sub style)
    # out -> a3,a4,a5,a6
    # borrow -> a1
    # ----------------------------

    # limb0
    sub     a3, t4, t0
    sltu    a1, t4, t0

    # limb1
    sub     t4, t5, t1           # diff
    sltu    t5, t5, t1           # borrow1
    sub     a4, t4, a1           # diff2
    sltu    t0, t4, a1           # borrow2   (reuse t0 as tmp)
    or      a1, t5, t0

    # limb2
    sub     t4, t6, t2
    sltu    t5, t6, t2
    sub     a5, t4, a1
    sltu    t0, t4, a1
    or      a1, t5, t0

    # limb3
    sub     t4, a7, t3
    sltu    t5, a7, t3
    sub     a6, t4, a1
    sltu    t0, t4, a1
    # or a1, t5, t0  # final borrow ignored

    # ----------------------------
    /*if a==0 => r=0 else r=t*/
    # ----------------------------
    and     a3, a3, a2
    and     a4, a4, a2
    and     a5, a5, a2
    and     a6, a6, a2

    sd      a3,   0(a0)
    sd      a4,   8(a0)
    sd      a5,  16(a0)
    sd      a6,  24(a0)

    ret
    .size fp_neg, .-fp_neg

    
