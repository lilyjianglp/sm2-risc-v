    .text
    .align  2
    .globl  fp_reduce
    .type   fp_reduce, @function

/* ============================================================
 * SM2 prime p:
 * p = FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF
 * limbs little-endian:
 *   p0 = 0xFFFFFFFFFFFFFFFF
 *   p1 = 0xFFFFFFFF00000000
 *   p2 = 0xFFFFFFFFFFFFFFFF
 *   p3 = 0xFFFFFFFEFFFFFFFF
 *
 * Reduce 512-bit T[8] -> r[4]
 * Using u384 accumulator x[0..5] in a2..a7:
 *   x0..x3 = T0..T3, x4=x5=0
 *
 * Fold relations:
 *   2^256 ≡ 1 + 2^224 + 2^96  - 2^64
 *   2^320 ≡ 1 + 2^32  + 2^160 + 2^224
 *   2^384 ≡ 1 + 2^32  + 2^96  + 2^128 + 2^225
 *   2^448 ≡ 2 + 2^33  + 2^96  + 2^129 + 2^160 + 2^192 + 2^225 - 2^64
 *
 * Tail rounds (fixed 2, unrolled):
 *   fold x4 as *2^256 and x5 as *2^320, repeat once
 *
 * Final: cond-sub p twice (constant-time)
 *
 * ABI:
 *   a0 = fp_t* r
 *   a1 = uint64_t* T
 *   uses a2..a7 as x[0..5]
 *   uses t0..t6 as scratch
 * ============================================================ */
 


/* ----- carry/borrow helpers (carry/borrow in t6) ----- */
/* carry-in/out in t6 (0/1) */
.macro ADDC64 rd, ra, rb
    add     \rd, \ra, \rb          /* sum */
    sltu    t4,  \rd, \ra          /* carry1 -> t4 */

    add     t5,  \rd, t6           /* sum2 = sum + cin */
    sltu    t6,  t5, \rd           /* carry2 -> t6 */

    mv      \rd, t5
    or      t6,  t4, t6            /* carry-out */
.endm

.macro ADD_CARRY_ONLY reg
    add     t5,  \reg, t6
    sltu    t6,  t5, \reg
    mv      \reg, t5
.endm

/* borrow-in/out in t6 (0/1) */
.macro SUBB64 rd, ra, rb
    sub     \rd, \ra, \rb          /* diff */
    sltu    t4,  \ra, \rb          /* borrow1 -> t4 */

    sub     t5,  \rd, t6           /* diff2 = diff - bin */
    sltu    t6,  \rd, t6           /* borrow2 -> t6 */

    mv      \rd, t5
    or      t6,  t4, t6            /* borrow-out */
.endm

.macro SUB_BORROW_ONLY reg
    sub     t5,  \reg, t6          /* reg2 = reg - bin */
    sltu    t6,  \reg, t6          /* borrow = (reg < bin) */
    mv      \reg, t5
.endm
/* ----- u384 add shifts: x[0..5] in a2..a7 ----- */

/* x += v<<0 */
.macro U384_ADD0 v
    li      t6, 0
    add     a2, a2, \v
    sltu    t6, a2, \v
    ADD_CARRY_ONLY a3
    ADD_CARRY_ONLY a4
    ADD_CARRY_ONLY a5
    ADD_CARRY_ONLY a6
    ADD_CARRY_ONLY a7
.endm

/* x += v<<32  -> (v<<32) into x0, (v>>32) into x1 */
.macro U384_ADD32 v
    slli    t2, \v, 32
    srli    t3, \v, 32
    li      t6, 0

    add     a2, a2, t2
    sltu    t6, a2, t2

    add     a3, a3, t3
    sltu    t0, a3, t3
    add     a3, a3, t6
    sltu    t6, a3, t6
    or      t6, t6, t0

    ADD_CARRY_ONLY a4
    ADD_CARRY_ONLY a5
    ADD_CARRY_ONLY a6
    ADD_CARRY_ONLY a7
.endm

/* x += v<<64 -> add v to x1 */
.macro U384_ADD64 v
    li      t6, 0
    add     a3, a3, \v
    sltu    t6, a3, \v
    ADD_CARRY_ONLY a4
    ADD_CARRY_ONLY a5
    ADD_CARRY_ONLY a6
    ADD_CARRY_ONLY a7
.endm

/* x += v<<96 -> (v<<32) into x1, (v>>32) into x2 */
.macro U384_ADD96 v
    slli    t2, \v, 32
    srli    t3, \v, 32
    li      t6, 0

    add     a3, a3, t2
    sltu    t6, a3, t2

    add     a4, a4, t3
    sltu    t0, a4, t3
    add     a4, a4, t6
    sltu    t6, a4, t6
    or      t6, t6, t0

    ADD_CARRY_ONLY a5
    ADD_CARRY_ONLY a6
    ADD_CARRY_ONLY a7
.endm

/* x += v<<128 -> add v to x2 */
.macro U384_ADD128 v
    li      t6, 0
    add     a4, a4, \v
    sltu    t6, a4, \v
    ADD_CARRY_ONLY a5
    ADD_CARRY_ONLY a6
    ADD_CARRY_ONLY a7
.endm

/* x += v<<129 -> (v<<1) into x2, (v>>63) into x3 */
.macro U384_ADD129 v
    slli    t2, \v, 1
    srli    t3, \v, 63
    li      t6, 0

    add     a4, a4, t2
    sltu    t6, a4, t2

    add     a5, a5, t3
    sltu    t0, a5, t3
    add     a5, a5, t6
    sltu    t6, a5, t6
    or      t6, t6, t0

    ADD_CARRY_ONLY a6
    ADD_CARRY_ONLY a7
.endm

/* x += v<<160 -> (v<<32) into x2, (v>>32) into x3 */
.macro U384_ADD160 v
    slli    t2, \v, 32
    srli    t3, \v, 32
    li      t6, 0

    add     a4, a4, t2
    sltu    t6, a4, t2

    add     a5, a5, t3
    sltu    t0, a5, t3
    add     a5, a5, t6
    sltu    t6, a5, t6
    or      t6, t6, t0

    ADD_CARRY_ONLY a6
    ADD_CARRY_ONLY a7
.endm

/* x += v<<192 -> add v to x3 */
.macro U384_ADD192 v
    li      t6, 0
    add     a5, a5, \v
    sltu    t6, a5, \v
    ADD_CARRY_ONLY a6
    ADD_CARRY_ONLY a7
.endm

/* x += v<<224 -> (v<<32) into x3, (v>>32) into x4 */
.macro U384_ADD224 v
    slli    t2, \v, 32
    srli    t3, \v, 32
    li      t6, 0

    add     a5, a5, t2
    sltu    t6, a5, t2

    add     a6, a6, t3
    sltu    t0, a6, t3
    add     a6, a6, t6
    sltu    t6, a6, t6
    or      t6, t6, t0

    ADD_CARRY_ONLY a7
.endm

/* x += v<<225 -> (v<<33) into x3, (v>>31) into x4 */
.macro U384_ADD225 v
    slli    t2, \v, 33
    srli    t3, \v, 31
    li      t6, 0

    add     a5, a5, t2
    sltu    t6, a5, t2

    add     a6, a6, t3
    sltu    t0, a6, t3
    add     a6, a6, t6
    sltu    t6, a6, t6
    or      t6, t6, t0

    ADD_CARRY_ONLY a7
.endm

/* x += v<<33 -> (v<<33) into x0, (v>>31) into x1 */
.macro U384_ADD33 v
    slli    t2, \v, 33
    srli    t3, \v, 31
    li      t6, 0

    add     a2, a2, t2
    sltu    t6, a2, t2

    add     a3, a3, t3
    sltu    t0, a3, t3
    add     a3, a3, t6
    sltu    t6, a3, t6
    or      t6, t6, t0

    ADD_CARRY_ONLY a4
    ADD_CARRY_ONLY a5
    ADD_CARRY_ONLY a6
    ADD_CARRY_ONLY a7
.endm


/* Real U384_SUB64: keep old a3 in t2 */
.macro U384_SUB64_REAL v
    li      t6, 0
    mv      t2, a3
    sub     a3, a3, \v
    sltu    t6, t2, \v        /* borrow */
    SUB_BORROW_ONLY a4
    SUB_BORROW_ONLY a5
    SUB_BORROW_ONLY a6
    SUB_BORROW_ONLY a7
.endm

/* ----- fold macros ----- */
.macro FOLD_2P256 k
    U384_ADD0    \k
    U384_ADD96   \k
    U384_ADD224  \k
    U384_SUB64_REAL \k
.endm

.macro FOLD_2P320 k
    U384_ADD0    \k
    U384_ADD32   \k
    U384_ADD160  \k
    U384_ADD224  \k
.endm

.macro FOLD_2P384 k
    U384_ADD0    \k
    U384_ADD32   \k
    U384_ADD96   \k
    U384_ADD128  \k
    U384_ADD225  \k
.endm

.macro FOLD_2P448 k
    /* +2*k */
    U384_ADD0    \k
    U384_ADD0    \k
    /* +2^33 */
    U384_ADD33   \k
    /* +2^96 */
    U384_ADD96   \k
    /* +2^129 */
    U384_ADD129  \k
    /* +2^160 */
    U384_ADD160  \k
    /* +2^192 */
    U384_ADD192  \k
    /* +2^225 */
    U384_ADD225  \k
    /* -2^64 */
    U384_SUB64_REAL \k
.endm

/* ----- cond-sub p (constant-time), applied to a2..a5 ----- */
.macro COND_SUB_P
    /* load p limbs (avoid t4/t5 because SUBB64 uses them as scratch now) */
    li      t2, -1                      /* p0 */
    li      t3, 0xFFFFFFFF00000000      /* p1 */
    li      a6, -1                      /* p2 */
    li      a7, 0xFFFFFFFEFFFFFFFF      /* p3 */

    /* tmp = x - p, borrow in t6 */
    li      t6, 0
    SUBB64  t0, a2, t2          /* tmp0 */
    SUBB64  t1, a3, t3          /* tmp1 */
    SUBB64  t2, a4, a6          /* tmp2 */
    SUBB64  t3, a5, a7          /* tmp3 */

    /* mask = -(borrow^1) */
    xori    t6, t6, 1
    neg     t6, t6
    not     t4, t6              /* ~mask (t4 is scratch OK) */

    /* select: x = (tmp & mask) | (x & ~mask) */
    and     t5, t0, t6
    and     a2, a2, t4
    or      a2, a2, t5

    and     t5, t1, t6
    and     a3, a3, t4
    or      a3, a3, t5

    and     t5, t2, t6
    and     a4, a4, t4
    or      a4, a4, t5

    and     t5, t3, t6
    and     a5, a5, t4
    or      a5, a5, t5
.endm


 
    
fp_reduce:
    addi    sp, sp, -32
    sd      s0,  0(sp)
    sd      s1,  8(sp)
    sd      s2, 16(sp)
    sd      s3, 24(sp)





    /* x0..x3 = T0..T3 */
    ld      a2,   0(a1)
    ld      a3,   8(a1)
    ld      a4,  16(a1)
    ld      a5,  24(a1)
    li      a6,  0
    li      a7,  0
    
    
 
    

    /* load high limbs t4..t7 */
    ld      s0,  32(a1)     /* T4 */
    ld      s1,  40(a1)     /* T5 */
    ld      s2,  48(a1)     /* T6 */
    ld      s3,  56(a1)     /* T7 */

    /* fold high limbs */
    FOLD_2P256 s0
 
    
    FOLD_2P320 s1
  
    
    FOLD_2P384 s2
  
    
    
    
    FOLD_2P448 s3
  

    /* tail round #1 */
    mv      s0, a6
    mv      s1, a7
    li      a6, 0
    li      a7, 0
    FOLD_2P256 s0
    FOLD_2P320 s1
   
    
    
    
    

    /* tail round #2 */
    mv      s0, a6
    mv      s1, a7
    li      a6, 0
    li      a7, 0
    FOLD_2P256 s0
    FOLD_2P320 s1
   
    
    
  
  
    /* final cond-sub p twice */
    COND_SUB_P
    COND_SUB_P
 

    /* store r->v[0..3] */
    sd      a2,  0(a0)
    sd      a3,  8(a0)
    sd      a4, 16(a0)
    sd      a5, 24(a0)
    
    ld      s0,  0(sp)
    ld      s1,  8(sp)
    ld      s2, 16(sp)
    ld      s3, 24(sp)
    addi    sp, sp, 32
    
    
    ret

    .size fp_reduce, .-fp_reduce
