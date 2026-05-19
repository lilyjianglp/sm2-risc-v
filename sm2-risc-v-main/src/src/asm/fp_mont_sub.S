    .text
    .align  2
    .globl  fp_mont_sub
    .type   fp_mont_sub, @function

/* ============================================================
 * fp_mont_sub(fp_t* r/a0, const fp_t* a/a1, const fp_t* b/a2)
 *
 * Constant-time Montgomery-field subtraction for 256-bit values:
 *   t = a - b
 *   r = borrow ? t + p : t
 *
 * Applied optimizations:
 *   - interleaved a[i]/b[i] loads, with constant build mixed in
 *   - subtract keeps cross-limb pipelining (OPT-C)
 *   - add-back precomputes cond-add on limbs 0/1/2/3
 *   - stores remain delayed off the critical carry path
 * ============================================================ */

fp_mont_sub:

    /* --------------------------------------------------------
     * interleaved loads with mixed constant build
     *   t6 = -1
     *   t4 = p[1] = 0xFFFFFFFF00000000
     *   t5 = p[3] = 0xFFFFFFFEFFFFFFFF = p[1] - 1
     * -------------------------------------------------------- */
    ld      a3,  0(a1)          /* a[0] */
    ld      t0,  0(a2)          /* b[0] */
    li      t6, -1
    ld      a4,  8(a1)          /* a[1] */
    ld      t1,  8(a2)          /* b[1] */
    slli    t4, t6, 32
    ld      a5, 16(a1)          /* a[2] */
    ld      t2, 16(a2)          /* b[2] */
    addi    t5, t4, -1
    ld      a6, 24(a1)          /* a[3] */
    ld      t3, 24(a2)          /* b[3] */

    /* --------------------------------------------------------
     * subtract stage with cross-limb pipelining
     * -------------------------------------------------------- */

    /* limb 0 */
    sltu    a7, a3, t0
    sub     a3, a3, t0

    /* precompute limb1 / limb2 partial subtracts */
    sltu    a1, a4, t1
    sub     a4, a4, t1
    sltu    t0, a5, t2
    sub     a5, a5, t2

    /* finalize limb1 */
    sltu    a2, a4, a7
    sub     a4, a4, a7
    or      a7, a1, a2

    /* precompute limb3 while finalizing limb2 */
    sltu    a1, a6, t3
    sub     a6, a6, t3
    sltu    a2, a5, a7
    sub     a5, a5, a7
    or      a7, t0, a2

    /* finalize limb3 */
    sltu    a2, a6, a7
    sub     a6, a6, a7
    or      a7, a1, a2          /* final borrow */

    /* --------------------------------------------------------
     * conditional add-back values
     *   t0 = mask = borrow ? -1 : 0
     *   t1 = cond_p1
     *   t3 = cond_p3
     * -------------------------------------------------------- */
    czero.eqz  t0, t6, a7
    and        t1, t4, t0
    and        t3, t5, t0

    /* --------------------------------------------------------
     * add-back scheduling
     *
     * First precompute cond-add on all limbs:
     *   a3 += t0
     *   a4 += t1
     *   a5 += t0
     *   a6 += t3
     *
     * Then run carry chain:
     *   carry0 from limb0 in a1
     *   carry(cond-add limb1) in a2
     *   carry(cond-add limb2) in t2
     *
     * This exposes more ILP than the strictly sequential form.
     * -------------------------------------------------------- */

    /* precompute cond-add on all limbs */
    add     a3, a3, t0
    add     a4, a4, t1
    add     a5, a5, t0
    add     a6, a6, t3

    /* carry from limb0 cond-add */
    sltu    a1, a3, t0

    /* precompute carry from limb1/2 cond-add */
    sltu    a2, a4, t1
    sltu    t2, a5, t0

    /* finalize limb1 */
    add     a4, a4, a1
    sd      a3,  0(a0)
    sltu    a1, a4, a1
    or      a1, a2, a1

    /* finalize limb2 */
    add     a5, a5, a1
    sd      a4,  8(a0)
    sltu    a1, a5, a1
    or      a1, t2, a1

    /* finalize limb3 */
    add     a6, a6, a1
    sd      a5, 16(a0)
    sd      a6, 24(a0)

    ret

.size fp_mont_sub, .-fp_mont_sub
