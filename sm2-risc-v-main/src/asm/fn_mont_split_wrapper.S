    .option norvc
    .text
    .align 2

    .globl fn_mont_mul_ct
    .type  fn_mont_mul_ct, @function
fn_mont_mul_ct:
    addi    sp, sp, -112
    sd      s2, 80(sp)
    sd      s1, 88(sp)
    sd      s0, 96(sp)
    sd      ra, 104(sp)

    mv      s0, a0              # save r
    mv      s1, a1              # save a
    mv      s2, a2              # save b

    # call fn_mul_4x4_u512_ct(t, a, b)
    addi    a0, sp, 0           # &t[0]
    mv      a1, s1
    mv      a2, s2
    call    fn_mul_4x4_u512_ct

    # t[8] = 0
    sd      zero, 64(sp)

    # call fn_mont_reduce_ct(r, t)
    mv      a0, s0
    addi    a1, sp, 0
    call    fn_mont_reduce_ct

    ld      s2, 80(sp)
    ld      s1, 88(sp)
    ld      s0, 96(sp)
    ld      ra, 104(sp)
    addi    sp, sp, 112
    ret

    .size fn_mont_mul_ct, .-fn_mont_mul_ct
