    .section .text
    .align  2
    .globl  fp_mont_neg
    .type   fp_mont_neg, @function

/*
 * void fp_mont_neg(fp_t *r, const fp_t *a)
 *
 * SM2 p limbs (little-endian):
 *   p0 = 0xFFFFFFFFFFFFFFFF
 *   p1 = 0xFFFFFFFF00000000
 *   p2 = 0xFFFFFFFFFFFFFFFF
 *   p3 = 0xFFFFFFFEFFFFFFFF
 *
 * r = (a == 0) ? 0 : (p - a)
 *
 * a0 = r
 * a1 = a
 */

fp_mont_neg:
    # load a[0..3]
    ld      t0, 0(a1)              # a0
    ld      t1, 8(a1)              # a1
    ld      t2, 16(a1)             # a2
    ld      t3, 24(a1)             # a3

    # nm = all-ones if a != 0 else 0
    or      t4, t0, t1
    or      t4, t4, t2
    or      t4, t4, t3
    snez    t4, t4
    sub     t4, x0, t4             # t4 = 0xffff...ffff if a!=0 else 0

    # --------------------------------------------------
    # limb0: p0 - a0 = 0xFFFFFFFFFFFFFFFF - a0 = ~a0
    # borrow0 = 0
    # --------------------------------------------------
    not     t5, t0                 # r0
    li      t6, 0                  # borrow = 0

    # --------------------------------------------------
    # limb1: p1 - a1 - borrow
    # p1 = 0xFFFFFFFF00000000
    # --------------------------------------------------
    li      a2, 0xffffffff00000000
    sub     a3, a2, t1             # d1
    sltu    a4, a2, t1             # b1
    sub     t1, a3, t6             # r1
    sltu    t6, a3, t6             # b2
    or      t6, a4, t6             # borrow1

    # --------------------------------------------------
    # limb2: p2 - a2 - borrow
    # p2 = 0xFFFFFFFFFFFFFFFF
    # --------------------------------------------------
    not     a3, t2                 # ~a2
    sub     t2, a3, t6             # r2
    sltu    t6, a3, t6             # borrow2

    # --------------------------------------------------
    # limb3: p3 - a3 - borrow
    # p3 = 0xFFFFFFFEFFFFFFFF
    # --------------------------------------------------
    li      a2, 0xfffffffeffffffff
    sub     a3, a2, t3             # d3
    sltu    a4, a2, t3             # b1
    sub     t3, a3, t6             # r3
    sltu    t6, a3, t6             # b2
    or      t6, a4, t6             # final borrow (unused)

   
    and     t5, t5, t4
    and     t1, t1, t4
    and     t2, t2, t4
    and     t3, t3, t4

    # store
    sd      t5, 0(a0)
    sd      t1, 8(a0)
    sd      t2, 16(a0)
    sd      t3, 24(a0)
    ret

    .size fp_mont_neg, .-fp_mont_neg
