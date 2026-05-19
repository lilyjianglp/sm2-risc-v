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


