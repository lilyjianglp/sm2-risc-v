    .text
    .align  2
    .globl  fn_mont_reduce_ct
    .type   fn_mont_reduce_ct, @function
fn_mont_reduce_ct:
    /* save s-registers (callee-saved) */
    addi    sp, sp, -48
    sd      s0, 0(sp)
    sd      s1, 8(sp)
    sd      s2, 16(sp)
    sd      s3, 24(sp)
    sd      s4, 32(sp)
    sd      ra, 40(sp)

    /* load constants */
    la      t6, .L_n
    ld      s0, 0(t6)              /* n0 */
    ld      s1, 8(t6)              /* n1 */
    ld      s2, 16(t6)             /* n2 */
    ld      s3, 24(t6)             /* n3 */

    la      t0, .L_n0prime
    ld      s4, 0(t0)              /* n0prime */

    la      t5, .L_mont_one        /* &mont_one[0] */

    /* for (i=0; i<4; i++) */
    li      t0, 0                  /* i */
.L_i_loop:
    slli    t1, t0, 3
    add     t2, a1, t1             /* t2 = &t[i] */

    /* m = (t[i] * n0prime) low64 */
    ld      a2, 0(t2)
    mul     a3, a2, s4             /* m */

    li      a4, 0                  /* carry = 0 */

    /* --------------------------------------------------------- */
    /* j=0 : make a4 the new high/carry directly, no mv a4,a7    */
    /* --------------------------------------------------------- */
    mul     a6, a3, s0             /* lo  = m*n0 */
    mulhu   a4, a3, s0             /* hi  = m*n0 high, directly into carry reg */
    ld      a2, 0(t2)
    add     a6, a6, a2
    sltu    a2, a6, a2
    add     a4, a4, a2
    sd      a6, 0(t2)

    /* --------------------------------------------------------- */
    /* j=1 */
    /* a4 enters as previous carry, exits as new carry           */
    /* --------------------------------------------------------- */
    mul     a6, a3, s1
    mulhu   a7, a3, s1
    ld      a2, 8(t2)
    add     a6, a6, a2
    sltu    a2, a6, a2
    add     a7, a7, a2
    add     a6, a6, a4
    sltu    a2, a6, a4
    add     a4, a7, a2             /* new carry directly to a4 */
    sd      a6, 8(t2)

    /* --------------------------------------------------------- */
    /* j=2 */
    /* --------------------------------------------------------- */
    mul     a6, a3, s2
    mulhu   a7, a3, s2
    ld      a2, 16(t2)
    add     a6, a6, a2
    sltu    a2, a6, a2
    add     a7, a7, a2
    add     a6, a6, a4
    sltu    a2, a6, a4
    add     a4, a7, a2             /* new carry directly to a4 */
    sd      a6, 16(t2)

    /* --------------------------------------------------------- */
    /* j=3 */
    /* --------------------------------------------------------- */
    mul     a6, a3, s3
    mulhu   a7, a3, s3
    ld      a2, 24(t2)
    add     a6, a6, a2
    sltu    a2, a6, a2
    add     a7, a7, a2
    add     a6, a6, a4
    sltu    a2, a6, a4
    add     a4, a7, a2             /* new carry directly to a4 */
    sd      a6, 24(t2)

    /* fixed propagation to t[i+4..8] */
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
    blt     t0, t1, .L_i_loop

    /* r = t[4..7] */
    ld      t0, 32(a1)
    sd      t0, 0(a0)
    ld      t0, 40(a1)
    sd      t0, 8(a0)
    ld      t0, 48(a1)
    sd      t0, 16(a0)
    ld      t0, 56(a1)
    sd      t0, 24(a0)

    /* t8 mask */
    ld      t0, 64(a1)
    andi    t0, t0, 1
    neg     t0, t0

    /* load r */
    ld      t1, 0(a0)
    ld      t2, 8(a0)
    ld      t3, 16(a0)
    ld      t4, 24(a0)

    /* add mont_one & mask */
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

    /* -------- cond_sub_n #1 using s0..s3 -------- */
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

    /* -------- cond_sub_n #2 (same) -------- */
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

    /* restore and return */
    ld      s0, 0(sp)
    ld      s1, 8(sp)
    ld      s2, 16(sp)
    ld      s3, 24(sp)
    ld      s4, 32(sp)
    ld      ra, 40(sp)
    addi    sp, sp, 48
    ret
    .size   fn_mont_reduce_ct, .-fn_mont_reduce_ct

    .section .rodata
    .align  3
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
