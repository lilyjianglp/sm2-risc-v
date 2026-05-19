    .text
    .align 2
    .globl pornin_inner_loop_31_asm
    .type  pornin_inner_loop_31_asm, @function

/*
 * void pornin_inner_loop_31_asm(
 *     uint64_t xa, uint64_t xb,
 *     int64_t *pf0, int64_t *pg0,
 *     int64_t *pf1, int64_t *pg1);
 *
 * RISC-V ABI:
 *   a0 = xa
 *   a1 = xb
 *   a2 = pf0
 *   a3 = pg0
 *   a4 = pf1
 *   a5 = pg1
 *
 * Output:
 *   *pf0 = f0
 *   *pg0 = g0
 *   *pf1 = f1
 *   *pg1 = g1
 *
 * Implements the same logic as pornin_inner_loop_31_c().
 */
pornin_inner_loop_31_asm:
    /* Save output pointers on stack; function is leaf. */
    addi    sp, sp, -32
    sd      a2, 0(sp)
    sd      a3, 8(sp)
    sd      a4, 16(sp)
    sd      a5, 24(sp)

    /*
     * State:
     *   a0 = xa
     *   a1 = xb
     *   t0 = f0
     *   t1 = g0
     *   t2 = f1
     *   t3 = g1
     *   t4 = loop counter
     */
    li      t0, 1          /* f0 = 1 */
    li      t1, 0          /* g0 = 0 */
    li      t2, 0          /* f1 = 0 */
    li      t3, 1          /* g1 = 1 */
    li      t4, 31

.L_inner31_loop:
    /*
     * odd = xa & 1
     * swap_mask = -(odd & (xa < xb))
     */
    andi    a2, a0, 1          /* a2 = odd */
    sltu    a3, a0, a1         /* a3 = xa < xb */
    sub     a5, x0, a2         /* a5 = odd_mask */
    and     a3, a3, a2         /* a3 = odd & (xa < xb) */
    sub     a4, x0, a3         /* a4 = swap_mask */

    /*
     * Conditional swaps with the same mask.
     *
     * Compute the independent xor diffs first, then mask/apply them.
     * This keeps the two ALU slots busier on in-order dual-issue cores.
     */
    xor     t5, a0, a1         /* diff_xa */
    xor     a6, t0, t2         /* diff_f */
    xor     a7, t1, t3         /* diff_g */
    and     t5, t5, a4
    and     a6, a6, a4
    and     a7, a7, a4
    xor     a0, a0, t5
    xor     a1, a1, t5
    xor     t0, t0, a6
    xor     t2, t2, a6
    xor     t1, t1, a7
    xor     t3, t3, a7

    /*
     * xa_even = xa >> 1
     * xa_odd  = (xa - xb) >> 1
     */
    srli    a6, a0, 1          /* a6 = xa_even */
    sub     a7, a0, a1
    srli    a7, a7, 1          /* a7 = xa_odd */

    /*
     * xa = odd ? xa_odd : xa_even
     * select with:
     *   xa = xa_even ^ ((xa_even ^ xa_odd) & odd_mask)
     */
    xor     t5, a6, a7
    and     t5, t5, a5
    xor     a0, a6, t5

    /*
     * f0_odd = f0 - f1
     * g0_odd = g0 - g1
     * f0 = odd ? f0_odd : f0
     * g0 = odd ? g0_odd : g0
     */
    sub     a6, t0, t2         /* a6 = f0_odd */
    sub     a7, t1, t3         /* a7 = g0_odd */

    xor     t5, t0, a6
    and     t5, t5, a5
    xor     t0, t0, t5

    xor     t5, t1, a7
    and     t5, t5, a5
    xor     t1, t1, t5

    /*
     * f1 *= 2
     * g1 *= 2
     */
    add     t2, t2, t2
    add     t3, t3, t3

    addi    t4, t4, -1
    bnez    t4, .L_inner31_loop

    /*
     * Store outputs.
     */
    ld      a2, 0(sp)
    ld      a3, 8(sp)
    ld      a4, 16(sp)
    ld      a5, 24(sp)

    sd      t0, 0(a2)
    sd      t1, 0(a3)
    sd      t2, 0(a4)
    sd      t3, 0(a5)

    addi    sp, sp, 32
    ret

    .size pornin_inner_loop_31_asm, .-pornin_inner_loop_31_asm
