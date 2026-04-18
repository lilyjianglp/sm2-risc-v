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