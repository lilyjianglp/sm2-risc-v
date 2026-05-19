    .text
    .align  2

    .globl  fp_mont_add
    .type   fp_mont_add, @function

/* ============================================================
 * fp_mont_add(fp_t* r/a0, const fp_t* a/a1, const fp_t* b/a2)
 *
 * Lazy experimental version.
 *
 * Input invariant for this fast path:
 *   a < p, b < p
 *
 * Output:
 *   r == a + b (mod p), but r is not always canonical.
 *
 * Instead of:
 *   t = a + b; if (carry || t >= p) t -= p
 *
 * this does:
 *   t = low_256(a + b)
 *   if carry: t += c, where c = 2^256 - p
 *
 * For SM2 p:
 *   c = [0000000000000001, 00000000FFFFFFFF,
 *        0000000000000000, 0000000100000000]
 *
 * This removes the full t-p borrow chain and the final select.
 * Do not feed repeatedly-lazy values here unless the caller tracks range
 * or normalizes before the next operation that requires canonical inputs.
 * ============================================================ */
fp_mont_add:
    /* interleaved loads */
    ld      a3,  0(a1)          /* a[0] */
    ld      t0,  0(a2)          /* b[0] */
    ld      a4,  8(a1)          /* a[1] */
    ld      t1,  8(a2)          /* b[1] */
    ld      a5, 16(a1)          /* a[2] */
    ld      t2, 16(a2)          /* b[2] */
    ld      a6, 24(a1)          /* a[3] */
    ld      t3, 24(a2)          /* b[3] */

    /* t = a + b
     * keep:
     *   t5 = t[0], t6 = t[1], a7 = t[2], a3 = t[3]
     *   t2 = final carry
     */

    /* limb0 */
    add     t5, a3, t0
    sltu    t0, t5, a3

    /* limb1 */
    add     t6, a4, t1
    sltu    t1, t6, a4
    add     t6, t6, t0
    sltu    t0, t6, t0
    or      t0, t1, t0

    /* limb2 */
    add     a7, a5, t2
    sltu    t1, a7, a5
    add     a7, a7, t0
    sltu    t0, a7, t0
    or      t0, t1, t0

    /* limb3 */
    add     a3, a6, t3
    sltu    t1, a3, a6
    add     a3, a3, t0
    sltu    t0, a3, t0
    or      t2, t1, t0          /* final carry: 0 or 1 */

    /* mask = -carry */
    neg     t2, t2

    /* t += c & mask
     * c0 = 1
     */
    andi    t0, t2, 1
    add     t5, t5, t0
    sltu    t4, t5, t0

    /* c1 = 0x00000000FFFFFFFF */
    li      t0, -1
    srli    t0, t0, 32
    and     t0, t0, t2
    add     t6, t6, t0
    sltu    t1, t6, t0
    add     t6, t6, t4
    sltu    t4, t6, t4
    or      t4, t4, t1

    /* c2 = 0 */
    add     a7, a7, t4
    sltu    t4, a7, t4

    /* c3 = 0x0000000100000000 */
    li      t0, 1
    slli    t0, t0, 32
    and     t0, t0, t2
    add     a3, a3, t0
    add     a3, a3, t4

    sd      t5,  0(a0)
    sd      t6,  8(a0)
    sd      a7, 16(a0)
    sd      a3, 24(a0)

    ret

.size fp_mont_add, .-fp_mont_add
