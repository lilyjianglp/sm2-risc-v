	.file	"fp_montgomery.c"
	.option nopic
	.attribute arch, "rv64i2p1_m2p0_a2p1_f2p2_d2p2_c2p0_zicond_zicsr2p0_zifencei2p0_zba1p0_zbb1p0"
	.attribute unaligned_access, 1
	.attribute stack_align, 16
	.text
	.align	1
	.type	fp_mont_reduce, @function
fp_mont_reduce:
.LFB24:
	.cfi_startproc
	addi	sp,sp,-32
	.cfi_def_cfa_offset 32
	li	t5,-1
	li	t3,-2
	sd	s0,24(sp)
	sd	s1,16(sp)
	addi	a6,a1,32
	sd	s2,8(sp)
	slli	t4,t5,32
	sd	s3,0(sp)
	.cfi_offset 8, -8
	.cfi_offset 9, -16
	.cfi_offset 18, -24
	.cfi_offset 19, -32
	rori	t3,t3,32
	addi	a7,a1,64
	li	t1,0
	li	t6,4
.L6:
	ld	a2,-32(a6)
	sd	zero,-32(a6)
	ld	a4,-24(a6)
	ld	a5,-16(a6)
	mul	t2,a2,t4
	neg	s1,a2
	mulhu	s3,a2,t4
	ld	a3,-8(a6)
	sub	a5,a5,a2
	mul	t0,a2,t3
	mulhu	s2,a2,t5
	sltu	s1,a5,s1
	add	a4,t2,a4
	mulhu	s0,a2,t3
	sltu	t2,a4,t2
	add	a2,a4,a2
	add	t2,t2,s3
	sltu	a4,a2,a4
	add	a4,a4,t2
	add	t2,s1,s2
	add	a4,a5,a4
	add	a3,t0,a3
	sltu	a5,a4,a5
	sltu	t0,a3,t0
	add	a5,a5,t2
	add	t0,t0,s0
	add	a5,a3,a5
	sd	a2,-24(a6)
	sltu	a3,a5,a3
	sd	a4,-16(a6)
	sd	a5,-8(a6)
	add	a3,a3,t0
	beq	a3,zero,.L2
	mv	a5,a6
	j	.L5
.L22:
	addi	a5,a5,8
	beq	a4,zero,.L2
.L5:
	ld	a2,0(a5)
	add	a4,a2,a3
	li	a3,1
	sd	a4,0(a5)
	sltu	a4,a4,a2
	bne	a7,a5,.L22
.L2:
	addi	t1,t1,1
	addi	a6,a6,8
	bne	t1,t6,.L6
	ld	t4,32(a1)
	li	a7,1
	slli	a7,a7,32
	li	t5,-1
	slli	t5,t5,32
	addi	a6,a7,1
	sd	t4,0(a0)
	li	a5,-2
	ld	a4,40(a1)
	rori	a5,a5,32
	sd	a4,8(a0)
	ld	t3,48(a1)
	sd	t3,16(a0)
	ld	t0,56(a1)
	sd	t0,24(a0)
	ld	a2,64(a1)
	andi	a2,a2,1
	add	s0,t4,a2
	neg	a2,a2
	add.uw	a1,a2,a4
	sltu	t4,s0,t4
	add	t4,a1,t4
	addi	t6,s0,1
	sltu	a3,a1,a4
	sltu	a4,t4,a1
	snez	s2,t6
	add	a1,t4,a7
	or	a3,a3,a4
	sltu	s1,t4,t5
	sgtu	a4,s2,a1
	andi	a3,a3,0xff
	add	a3,t3,a3
	or	s1,s1,a4
	addi	t1,a3,1
	andi	s1,s1,0xff
	and	a2,a2,a7
	sgtu	t2,s1,t1
	sltu	t3,a3,t3
	add	a2,a2,t0
	snez	t0,t1
	add	a4,a2,t3
	or	t0,t0,t2
	add	t2,a4,a6
	andi	t0,t0,0xff
	sltu	a2,a4,a5
	sgtu	s3,t0,t2
	sub	t1,t1,s1
	or	a2,a2,s3
	sub	t3,a1,s2
	xori	a2,a2,1
	sub	t2,t2,t0
	andi	a2,a2,0xff
	ld	s2,8(sp)
	.cfi_restore 18
	neg	s1,a2
	addi	a2,a2,-1
	and	s0,s0,a2
	and	t6,t6,s1
	or	t6,t6,s0
	and	t4,t4,a2
	and	t3,t3,s1
	addi	t0,t6,1
	or	t3,t3,t4
	snez	s0,t0
	add	a7,t3,a7
	and	a3,a3,a2
	sgtu	t4,s0,a7
	sltu	a1,t3,t5
	and	t1,t1,s1
	or	a1,a1,t4
	or	t1,t1,a3
	andi	a1,a1,0xff
	addi	a3,t1,1
	and	a4,a4,a2
	sgtu	t5,a1,a3
	and	t4,t2,s1
	snez	a2,a3
	or	t4,t4,a4
	or	a2,a2,t5
	add	a4,t4,a6
	andi	a2,a2,0xff
	sltu	a5,t4,a5
	sgtu	a6,a2,a4
	sub	a7,a7,s0
	or	a5,a5,a6
	sub	a3,a3,a1
	xori	a5,a5,1
	sub	a4,a4,a2
	andi	a5,a5,0xff
	ld	s0,24(sp)
	.cfi_restore 8
	neg	a1,a5
	addi	a5,a5,-1
	and	t6,a5,t6
	and	t3,a5,t3
	and	t1,a5,t1
	and	t0,t0,a1
	and	a2,a7,a1
	and	a3,a3,a1
	and	a4,a4,a1
	and	a5,a5,t4
	or	t0,t0,t6
	or	a2,a2,t3
	or	a3,a3,t1
	or	a5,a4,a5
	sd	t0,0(a0)
	ld	s1,16(sp)
	.cfi_restore 9
	sd	a2,8(a0)
	ld	s3,0(sp)
	.cfi_restore 19
	sd	a3,16(a0)
	sd	a5,24(a0)
	addi	sp,sp,32
	.cfi_def_cfa_offset 0
	jr	ra
	.cfi_endproc
.LFE24:
	.size	fp_mont_reduce, .-fp_mont_reduce
	.align	1
	.type	comba_sqr, @function
comba_sqr:
.LFB23:
	.cfi_startproc
	ld	t4,8(a1)
	addi	sp,sp,-32
	.cfi_def_cfa_offset 32
	ld	a6,0(a1)
	sd	s0,24(sp)
	ld	a4,16(a1)
	sd	s2,8(sp)
	.cfi_offset 8, -8
	.cfi_offset 18, -24
	mul	s0,t4,t4
	sd	s3,0(sp)
	mul	t1,t4,a6
	ld	a1,24(a1)
	mulhu	a5,t4,a6
	sd	s1,16(sp)
	.cfi_offset 19, -32
	.cfi_offset 9, -16
	mulhu	a7,a6,a6
	mul	a3,a4,a6
	mulhu	t6,a4,a6
	srli	t3,t1,63
	mulhu	s2,t4,t4
	sh1add	a2,t1,a7
	slli	t1,a5,1
	or	t1,t1,t3
	sltu	a7,a2,a7
	add	a7,t1,a7
	srli	t3,a3,63
	sltu	t1,a7,t1
	slli	s3,t6,1
	srli	a5,a5,63
	or	s3,s3,t3
	sh1add	a7,a3,a7
	add	a5,a5,t1
	slli	a3,a3,1
	add	a5,s3,a5
	sltu	t5,a7,a3
	add	a7,s0,a7
	add	t2,t5,a5
	mulhu	t0,a1,a6
	mul	t1,a1,a6
	add	t3,s2,t2
	sltu	s0,a7,s0
	sltu	a3,a5,s3
	add	s1,s0,t3
	sltu	a5,t2,t5
	sltu	s0,s1,s0
	mul	t5,a4,t4
	mulhu	t2,a4,t4
	or	a3,a3,a5
	sltu	t3,t3,s2
	srli	t6,t6,63
	or	t3,t3,s0
	andi	a3,a3,0xff
	add	a3,a3,t6
	slli	a5,t0,1
	srli	t6,t1,63
	andi	t3,t3,0xff
	add	t3,t3,a3
	sh1add	t1,t1,s1
	or	a5,a5,t6
	sltu	s1,t1,s1
	add	a5,t3,a5
	srli	a3,t5,63
	slli	s2,t2,1
	add	s3,s1,a5
	or	s2,s2,a3
	sh1add	t1,t5,t1
	slli	t5,t5,1
	mul	s0,a1,t4
	mulhu	a3,a1,t4
	add	t6,s2,s3
	sltu	t5,t1,t5
	sltu	s3,s3,s1
	add	t4,t5,t6
	sltu	a5,a5,t3
	sltu	t5,t4,t5
	srli	t3,t0,63
	or	a5,a5,s3
	sltu	t6,t6,s2
	srli	t2,t2,63
	mul	t0,a4,a4
	add	t2,t2,t3
	or	t6,t6,t5
	andi	a5,a5,0xff
	mulhu	s1,a4,a4
	srli	t3,s0,63
	add	a5,a5,t2
	andi	t6,t6,0xff
	slli	s2,a3,1
	or	s2,s2,t3
	sh1add	t4,s0,t4
	add	a5,a5,t6
	slli	t5,s0,1
	sltu	t5,t4,t5
	add	a5,s2,a5
	add	s0,t5,a5
	add	t4,t0,t4
	mul	t2,a1,a4
	add	t3,s1,s0
	mulhu	a4,a1,a4
	sltu	t0,t4,t0
	add	t6,t0,t3
	sltu	s0,s0,t5
	sltu	a5,a5,s2
	sltu	t3,t3,s1
	or	a5,a5,s0
	sltu	t0,t6,t0
	srli	a3,a3,63
	or	t3,t3,t0
	andi	a5,a5,0xff
	srli	t5,t2,63
	add	a5,a5,a3
	andi	t3,t3,0xff
	slli	a3,a4,1
	add	t3,t3,a5
	sh1add	t2,t2,t6
	or	a5,a3,t5
	mul	t0,a1,a1
	add	a5,t3,a5
	sltu	a3,t2,t6
	mulhu	a1,a1,a1
	add	t5,a3,a5
	sltu	a5,a5,t3
	sltu	a3,t5,a3
	mul	a6,a6,a6
	or	a5,a5,a3
	srli	a4,a4,63
	add	a3,t0,t5
	add	a4,a4,a1
	andi	a5,a5,0xff
	ld	s0,24(sp)
	.cfi_restore 8
	add	a5,a5,a4
	sltu	t0,a3,t0
	add	a5,a5,t0
	sd	a6,0(a0)
	sd	a2,8(a0)
	ld	s1,16(sp)
	.cfi_restore 9
	sd	a7,16(a0)
	ld	s2,8(sp)
	.cfi_restore 18
	sd	t1,24(a0)
	ld	s3,0(sp)
	.cfi_restore 19
	sd	t4,32(a0)
	sd	t2,40(a0)
	sd	a3,48(a0)
	sd	a5,56(a0)
	addi	sp,sp,32
	.cfi_def_cfa_offset 0
	jr	ra
	.cfi_endproc
.LFE23:
	.size	comba_sqr, .-comba_sqr
	.align	1
	.globl	fp_cmp
	.type	fp_cmp, @function
fp_cmp:
.LFB7:
	.cfi_startproc
	li	a4,24
	li	a6,0
	li	a7,0
	li	t3,-8
.L74:
	add	a5,a0,a4
	add	a2,a1,a4
	ld	a3,0(a5)
	or	a5,a7,a6
	ld	t1,0(a2)
	xori	a5,a5,1
	andi	a5,a5,1
	addi	a4,a4,-8
	sgtu	a2,a3,t1
	sltu	a3,a3,t1
	and	a2,a5,a2
	and	a5,a5,a3
	or	a7,a7,a2
	or	a6,a6,a5
	bne	a4,t3,.L74
	subw	a0,a7,a6
	ret
	.cfi_endproc
.LFE7:
	.size	fp_cmp, .-fp_cmp
	.align	1
	.globl	fp_copy
	.type	fp_copy, @function
fp_copy:
.LFB8:
	.cfi_startproc
	ld	a2,0(a1)
	ld	a3,8(a1)
	ld	a4,16(a1)
	ld	a5,24(a1)
	sd	a2,0(a0)
	sd	a3,8(a0)
	sd	a4,16(a0)
	sd	a5,24(a0)
	ret
	.cfi_endproc
.LFE8:
	.size	fp_copy, .-fp_copy
	.align	1
	.globl	fp_set_zero
	.type	fp_set_zero, @function
fp_set_zero:
.LFB9:
	.cfi_startproc
	sd	zero,0(a0)
	sd	zero,8(a0)
	sd	zero,16(a0)
	sd	zero,24(a0)
	ret
	.cfi_endproc
.LFE9:
	.size	fp_set_zero, .-fp_set_zero
	.align	1
	.globl	fp_set_one
	.type	fp_set_one, @function
fp_set_one:
.LFB10:
	.cfi_startproc
	li	a5,1
	sd	zero,8(a0)
	sd	a5,0(a0)
	sd	zero,16(a0)
	sd	zero,24(a0)
	ret
	.cfi_endproc
.LFE10:
	.size	fp_set_one, .-fp_set_one
	.align	1
	.globl	fp_is_zero
	.type	fp_is_zero, @function
fp_is_zero:
.LFB11:
	.cfi_startproc
	ld	a5,0(a0)
	ld	a2,8(a0)
	ld	a3,16(a0)
	ld	a4,24(a0)
	or	a0,a5,a2
	or	a0,a0,a3
	or	a0,a0,a4
	seqz	a0,a0
	ret
	.cfi_endproc
.LFE11:
	.size	fp_is_zero, .-fp_is_zero
	.align	1
	.globl	fp_is_equal
	.type	fp_is_equal, @function
fp_is_equal:
.LFB12:
	.cfi_startproc
	ld	a4,0(a1)
	ld	a5,0(a0)
	ld	a2,8(a0)
	ld	a7,8(a1)
	ld	a3,16(a0)
	xor	a5,a5,a4
	ld	a6,16(a1)
	ld	a4,24(a0)
	xor	a2,a2,a7
	ld	a1,24(a1)
	or	a0,a5,a2
	xor	a5,a3,a6
	or	a0,a0,a5
	xor	a5,a4,a1
	or	a0,a0,a5
	seqz	a0,a0
	ret
	.cfi_endproc
.LFE12:
	.size	fp_is_equal, .-fp_is_equal
	.align	1
	.globl	fp_from_bytes
	.type	fp_from_bytes, @function
fp_from_bytes:
.LFB13:
	.cfi_startproc
	addi	sp,sp,-64
	.cfi_def_cfa_offset 64
	lui	a7,%hi(__stack_chk_guard)
	sd	ra,56(sp)
	.cfi_offset 1, -8
	ld	a5, %lo(__stack_chk_guard)(a7)
	sd	a5, 40(sp)
	li	a5, 0
	addi	a6,sp,8
	addi	a1,a1,8
	li	a2,0
	li	t1,32
.L82:
	addi	a4,a1,-8
	li	a5,0
.L83:
	lbu	a3,0(a4)
	slli	a5,a5,8
	addi	a4,a4,1
	or	a5,a3,a5
	bne	a4,a1,.L83
	sd	a5,24(a6)
	addiw	a2,a2,8
	addi	a6,a6,-8
	addi	a1,a4,8
	bne	a2,t1,.L82
	ld	t5,8(sp)
	li	a5,1
	ld	t4,16(sp)
	slli	a5,a5,32
	ld	t3,24(sp)
	li	a1,-1
	addi	a6,t5,1
	slli	a1,a1,32
	snez	t0,a6
	add	a3,t4,a5
	sltu	a5,a3,t0
	sltu	a1,t4,a1
	or	a1,a1,a5
	ld	t6,32(sp)
	addi	a2,t3,1
	andi	a1,a1,0xff
	sltu	a4,a2,a1
	snez	t1,a2
	li	a5,1
	slli	a5,a5,32
	or	t1,t1,a4
	addi	a5,a5,1
	add	a4,t6,a5
	andi	t1,t1,0xff
	li	a5,-2
	rori	a5,a5,32
	sltu	t2,a4,t1
	sltu	a5,t6,a5
	or	a5,a5,t2
	sub	a2,a2,a1
	xori	a5,a5,1
	sub	a4,a4,t1
	andi	a5,a5,0xff
	sub	a3,a3,t0
	neg	t1,a5
	addi	a5,a5,-1
	and	a1,a6,t1
	and	t5,t5,a5
	and	t4,t4,a5
	and	t3,t3,a5
	and	a4,a4,t1
	and	a3,a3,t1
	and	a2,a2,t1
	and	t6,t6,a5
	or	a5,a4,t6
	or	a1,a1,t5
	or	a3,a3,t4
	or	a2,a2,t3
	ld	a6, 40(sp)
	ld	a4, %lo(__stack_chk_guard)(a7)
	xor	a4, a6, a4
	li	a6, 0
	sd	a1,8(sp)
	sd	a3,16(sp)
	sd	a2,24(sp)
	sd	a5,32(sp)
	sd	a1,0(a0)
	sd	a3,8(a0)
	sd	a2,16(a0)
	sd	a5,24(a0)
	bne	a4,zero,.L88
	ld	ra,56(sp)
	.cfi_remember_state
	.cfi_restore 1
	addi	sp,sp,64
	.cfi_def_cfa_offset 0
	jr	ra
.L88:
	.cfi_restore_state
	call	__stack_chk_fail
	.cfi_endproc
.LFE13:
	.size	fp_from_bytes, .-fp_from_bytes
	.align	1
	.globl	fp_to_bytes
	.type	fp_to_bytes, @function
fp_to_bytes:
.LFB14:
	.cfi_startproc
	li	a7,0
	li	a6,-8
	li	t1,-32
.L91:
	add	a5,a1,a7
	sub	a4,a0,a7
	ld	a2,24(a5)
	li	a5,56
.L90:
	srl	a3,a2,a5
	addiw	a5,a5,-8
	sb	a3,0(a4)
	addi	a4,a4,1
	bne	a5,a6,.L90
	addi	a7,a7,-8
	bne	a7,t1,.L91
	ret
	.cfi_endproc
.LFE14:
	.size	fp_to_bytes, .-fp_to_bytes
	.align	1
	.globl	fp_add
	.type	fp_add, @function
fp_add:
.LFB15:
	.cfi_startproc
	ld	t1,0(a1)
	li	a7,1
	ld	t3,0(a2)
	slli	a7,a7,32
	ld	a6,8(a1)
	addi	sp,sp,-16
	.cfi_def_cfa_offset 16
	ld	a3,8(a2)
	sd	s0,8(sp)
	.cfi_offset 8, -8
	add	t3,t1,t3
	ld	a4,16(a1)
	sltu	t1,t3,t1
	ld	a5,16(a2)
	add	a3,a6,a3
	addi	t4,t3,1
	add	t1,a3,t1
	sltu	a6,a3,a6
	sltu	a3,t1,a3
	ld	t6,24(a1)
	or	a6,a6,a3
	ld	a1,24(a2)
	add	a7,t1,a7
	add	a5,a4,a5
	snez	t5,t4
	andi	a6,a6,0xff
	li	a2,-1
	slli	a2,a2,32
	add	a6,a5,a6
	sltu	a2,t1,a2
	sltu	a3,a7,t5
	sltu	a4,a5,a4
	or	a3,a2,a3
	sltu	a5,a6,a5
	addi	t0,a6,1
	andi	a2,a3,0xff
	or	a4,a4,a5
	add	a1,t6,a1
	sltu	a3,t0,a2
	andi	a4,a4,0xff
	snez	t2,t0
	li	a5,1
	slli	a5,a5,32
	add	a4,a1,a4
	or	t2,t2,a3
	addi	a5,a5,1
	add	a3,a4,a5
	andi	t2,t2,0xff
	li	a5,-2
	rori	a5,a5,32
	sltu	s0,a3,t2
	sltu	a5,a4,a5
	or	a5,a5,s0
	sltu	t6,a1,t6
	sltu	a1,a4,a1
	xori	a5,a5,1
	or	t6,t6,a1
	sub	a2,t0,a2
	or	a5,a5,t6
	sub	a1,a7,t5
	andi	a5,a5,0xff
	sub	a3,a3,t2
	neg	a7,a5
	addi	a5,a5,-1
	ld	s0,8(sp)
	.cfi_restore 8
	and	t4,t4,a7
	and	t3,t3,a5
	and	a1,a1,a7
	and	t1,t1,a5
	and	a2,a2,a7
	and	a6,a6,a5
	and	a3,a3,a7
	and	a4,a4,a5
	or	t4,t4,t3
	or	a1,a1,t1
	or	a2,a2,a6
	or	a3,a3,a4
	sd	t4,0(a0)
	sd	a1,8(a0)
	sd	a2,16(a0)
	sd	a3,24(a0)
	addi	sp,sp,16
	.cfi_def_cfa_offset 0
	jr	ra
	.cfi_endproc
.LFE15:
	.size	fp_add, .-fp_add
	.align	1
	.globl	fp_sub
	.type	fp_sub, @function
fp_sub:
.LFB16:
	.cfi_startproc
	ld	a4,0(a2)
	li	t4,-2
	ld	a7,8(a2)
	rori	t4,t4,32
	ld	t1,0(a1)
	ld	a5,8(a1)
	ld	a3,16(a1)
	sltu	a6,t1,a4
	sub	t1,t1,a4
	sub	a4,a5,a7
	ld	t2,16(a2)
	sub	t3,a4,a6
	sltu	a4,a4,a6
	li	a6,-1
	slli	a6,a6,32
	add	a6,t3,a6
	sltu	a5,a5,a7
	snez	a7,t1
	or	a4,a5,a4
	add	a7,a6,a7
	ld	a5,24(a1)
	ld	t0,24(a2)
	andi	a1,a4,0xff
	sub	t6,a3,t2
	sltu	a4,a6,t3
	sltu	a6,a7,a6
	sub	t5,t6,a1
	sltu	a2,a3,t2
	sltu	t6,t6,a1
	or	a4,a4,a6
	addi	a1,t5,-1
	or	a2,a2,t6
	andi	a4,a4,0xff
	sub	a6,a5,t0
	andi	a2,a2,0xff
	add	a4,a1,a4
	sltu	t6,a6,a2
	sltu	a1,a4,a1
	snez	a3,t5
	sltu	a5,a5,t0
	sub	a2,a6,a2
	or	a3,a3,a1
	or	a5,a5,t6
	andi	a5,a5,0xff
	add	t4,a2,t4
	andi	a3,a3,0xff
	neg	a6,a5
	addi	a1,t1,-1
	addi	a5,a5,-1
	add	a3,a3,t4
	and	a1,a1,a6
	and	t1,t1,a5
	and	a7,a7,a6
	and	t3,t3,a5
	and	a4,a4,a6
	and	t5,t5,a5
	and	a3,a3,a6
	and	a2,a2,a5
	or	a1,a1,t1
	or	a7,a7,t3
	or	a4,a4,t5
	or	a3,a3,a2
	sd	a1,0(a0)
	sd	a7,8(a0)
	sd	a4,16(a0)
	sd	a3,24(a0)
	ret
	.cfi_endproc
.LFE16:
	.size	fp_sub, .-fp_sub
	.align	1
	.globl	fp_mont_set_zero
	.type	fp_mont_set_zero, @function
fp_mont_set_zero:
.LFB33:
	.cfi_startproc
	sd	zero,0(a0)
	sd	zero,8(a0)
	sd	zero,16(a0)
	sd	zero,24(a0)
	ret
	.cfi_endproc
.LFE33:
	.size	fp_mont_set_zero, .-fp_mont_set_zero
	.align	1
	.globl	fp_mont_set_one
	.type	fp_mont_set_one, @function
fp_mont_set_one:
.LFB18:
	.cfi_startproc
	lui	a5,%hi(.LANCHOR0)
	addi	a5,a5,%lo(.LANCHOR0)
	ld	a2,0(a5)
	ld	a3,8(a5)
	ld	a4,16(a5)
	ld	a5,24(a5)
	sd	a2,0(a0)
	sd	a3,8(a0)
	sd	a4,16(a0)
	sd	a5,24(a0)
	ret
	.cfi_endproc
.LFE18:
	.size	fp_mont_set_one, .-fp_mont_set_one
	.align	1
	.globl	fp_mont_add
	.type	fp_mont_add, @function
fp_mont_add:
.LFB35:
	.cfi_startproc
	tail	fp_add
	.cfi_endproc
.LFE35:
	.size	fp_mont_add, .-fp_mont_add
	.align	1
	.globl	fp_mont_sub
	.type	fp_mont_sub, @function
fp_mont_sub:
.LFB37:
	.cfi_startproc
	tail	fp_sub
	.cfi_endproc
.LFE37:
	.size	fp_mont_sub, .-fp_mont_sub
	.align	1
	.globl	fp_mont_neg
	.type	fp_mont_neg, @function
fp_mont_neg:
.LFB21:
	.cfi_startproc
	ld	a7,8(a1)
	li	a5,-1
	ld	a6,0(a1)
	slli	a5,a5,32
	ld	a2,16(a1)
	addi	a5,a5,1
	ld	a1,24(a1)
	sltu	a4,a7,a5
	or	a5,a6,a7
	xori	a4,a4,1
	or	a5,a5,a2
	li	a3,-2
	not	a2,a2
	or	a5,a5,a1
	rori	a3,a3,32
	sltu	t1,a2,a4
	sub	a3,a3,a1
	seqz	a5,a5
	li	a1,-1
	slli	a1,a1,32
	addi	a5,a5,-1
	sub	a1,a1,a7
	sub	a2,a2,a4
	sub	a3,a3,t1
	andn	a6,a5,a6
	and	a1,a1,a5
	and	a2,a2,a5
	and	a3,a3,a5
	sd	a6,0(a0)
	sd	a1,8(a0)
	sd	a2,16(a0)
	sd	a3,24(a0)
	ret
	.cfi_endproc
.LFE21:
	.size	fp_mont_neg, .-fp_mont_neg
	.align	1
	.globl	fp_mul_4x4_u512_ct
	.type	fp_mul_4x4_u512_ct, @function
fp_mul_4x4_u512_ct:
.LFB22:
	.cfi_startproc
	ld	a7,0(a2)
	addi	sp,sp,-96
	.cfi_def_cfa_offset 96
	ld	t0,8(a2)
	sd	s1,80(sp)
	ld	a6,0(a1)
	sd	s2,72(sp)
	ld	t2,8(a1)
	sd	s5,48(sp)
	ld	t6,16(a1)
	sd	s0,88(sp)
	mulhu	a5,a6,a7
	ld	t3,24(a1)
	mul	a4,t0,a6
	ld	t4,24(a2)
	mul	a3,t2,a7
	ld	t5,16(a2)
	mulhu	t1,t0,a6
	sd	s3,64(sp)
	mulhu	a1,t2,a7
	sd	s4,56(sp)
	add	a2,a4,a5
	.cfi_offset 9, -16
	.cfi_offset 18, -24
	.cfi_offset 21, -48
	.cfi_offset 8, -8
	.cfi_offset 19, -32
	.cfi_offset 20, -40
	mul	s1,t5,a6
	sltu	a4,a2,a4
	add	a2,a3,a2
	add	a4,a4,t1
	sltu	a3,a2,a3
	add	a4,a1,a4
	mulhu	s2,t5,a6
	add	t1,a3,a4
	mul	s5,t2,t0
	sltu	a5,a4,a1
	sltu	a3,t1,a3
	mulhu	s0,t2,t0
	or	a5,a5,a3
	mul	s3,t6,a7
	add	t1,s1,t1
	andi	a5,a5,0xff
	mulhu	s4,t6,a7
	sltu	s1,t1,s1
	add	a5,a5,s2
	add	a5,s1,a5
	add	t1,s5,t1
	add	a1,s0,a5
	sltu	s5,t1,s5
	sd	s6,40(sp)
	add	t1,s3,t1
	.cfi_offset 22, -56
	add	s6,s5,a1
	sltu	s3,t1,s3
	add	a4,s4,s6
	mul	s2,t4,a6
	add	a3,s3,a4
	sltu	a1,a1,s0
	sltu	s0,s6,s5
	mulhu	s6,t4,a6
	or	s0,a1,s0
	sltu	a1,a4,s4
	sltu	s4,a3,s3
	andi	a4,s0,0xff
	mul	s3,t5,t2
	sltu	a5,a5,s1
	or	a1,a1,s4
	sd	s7,32(sp)
	add	a5,a4,a5
	.cfi_offset 23, -64
	mulhu	s7,t5,t2
	andi	a1,a1,0xff
	sd	s9,16(sp)
	add	a3,s2,a3
	.cfi_offset 25, -80
	mul	s9,t6,t0
	add	a1,a1,a5
	sd	s8,24(sp)
	sltu	s2,a3,s2
	.cfi_offset 24, -72
	mulhu	s8,t6,t0
	add	a1,s6,a1
	mul	a5,t3,a7
	add	s4,s2,a1
	add	a3,s3,a3
	mulhu	a4,t3,a7
	add	s5,s7,s4
	sltu	s3,a3,s3
	add	a3,s9,a3
	add	s1,s3,s5
	sltu	s9,a3,s9
	add	s0,s8,s1
	sd	s10,8(sp)
	add	a3,a5,a3
	.cfi_offset 26, -88
	add	s10,s9,s0
	sd	s11,0(sp)
	.cfi_offset 27, -96
	sltu	s5,s5,s7
	sltu	s11,a3,a5
	sltu	s4,s4,s2
	add	a5,a4,s10
	sltu	s1,s1,s3
	sltu	a1,a1,s6
	mul	s2,t4,t2
	add	s3,s11,a5
	sltu	s10,s10,s9
	or	s1,s5,s1
	or	a1,a1,s4
	sltu	s0,s0,s8
	mulhu	t2,t4,t2
	andi	a1,a1,0xff
	or	s0,s0,s10
	sltu	a4,a5,a4
	andi	s1,s1,0xff
	sltu	a5,s3,s11
	mul	s6,t6,t5
	add	s1,s1,a1
	or	a4,a4,a5
	andi	s0,s0,0xff
	mulhu	s4,t6,t5
	add	s0,s0,s1
	andi	a4,a4,0xff
	mul	s1,t3,t0
	add	a1,s2,s3
	add	a4,a4,s0
	mulhu	s5,t3,t0
	sltu	s2,a1,s2
	add	a4,t2,a4
	add	s7,s2,a4
	add	a1,s6,a1
	add	t0,s4,s7
	sltu	s6,a1,s6
	add	s8,s6,t0
	add	a1,s1,a1
	add	a5,s5,s8
	sltu	s3,a1,s1
	mul	s0,t4,t6
	add	s1,s3,a5
	sltu	s8,s8,s6
	sltu	s7,s7,s2
	sltu	t0,t0,s4
	sltu	a4,a4,t2
	mulhu	t6,t4,t6
	sltu	s3,s1,s3
	or	t0,t0,s8
	or	a4,a4,s7
	sltu	a5,a5,s5
	mul	t2,t3,t5
	andi	a4,a4,0xff
	or	a5,a5,s3
	andi	t0,t0,0xff
	mulhu	t5,t3,t5
	add	t0,t0,a4
	andi	a5,a5,0xff
	add	s1,s0,s1
	add	a5,a5,t0
	sltu	s0,s1,s0
	add	a5,t6,a5
	add	s2,s0,a5
	add	t0,t2,s1
	add	a4,t5,s2
	mul	s1,t3,t4
	sltu	t2,t0,t2
	mulhu	t3,t3,t4
	sltu	s2,s2,s0
	add	t4,t2,a4
	sltu	a5,a5,t6
	sltu	a4,a4,t5
	or	a5,a5,s2
	sltu	t2,t4,t2
	mul	a6,a6,a7
	andi	a5,a5,0xff
	or	a4,a4,t2
	ld	s0,88(sp)
	.cfi_restore 8
	add	a7,s1,t4
	add	a5,a5,t3
	andi	a4,a4,0xff
	sltu	s1,a7,s1
	add	a5,a4,a5
	ld	s2,72(sp)
	.cfi_restore 18
	add	a5,a5,s1
	ld	s3,64(sp)
	.cfi_restore 19
	ld	s1,80(sp)
	.cfi_restore 9
	sd	a6,0(a0)
	sd	a2,8(a0)
	sd	t1,16(a0)
	sd	a3,24(a0)
	sd	a1,32(a0)
	sd	t0,40(a0)
	sd	a7,48(a0)
	sd	a5,56(a0)
	ld	s4,56(sp)
	.cfi_restore 20
	ld	s5,48(sp)
	.cfi_restore 21
	ld	s6,40(sp)
	.cfi_restore 22
	ld	s7,32(sp)
	.cfi_restore 23
	ld	s8,24(sp)
	.cfi_restore 24
	ld	s9,16(sp)
	.cfi_restore 25
	ld	s10,8(sp)
	.cfi_restore 26
	ld	s11,0(sp)
	.cfi_restore 27
	addi	sp,sp,96
	.cfi_def_cfa_offset 0
	jr	ra
	.cfi_endproc
.LFE22:
	.size	fp_mul_4x4_u512_ct, .-fp_mul_4x4_u512_ct
	.align	1
	.globl	fp_mont_reduce_ct
	.type	fp_mont_reduce_ct, @function
fp_mont_reduce_ct:
.LFB25:
	.cfi_startproc
	tail	fp_mont_reduce
	.cfi_endproc
.LFE25:
	.size	fp_mont_reduce_ct, .-fp_mont_reduce_ct
	.align	1
	.globl	fp_mont_mul_ct
	.type	fp_mont_mul_ct, @function
fp_mont_mul_ct:
.LFB26:
	.cfi_startproc
	addi	sp,sp,-176
	.cfi_def_cfa_offset 176
	sd	s0,160(sp)
	.cfi_offset 8, -16
	mv	s0,a0
	sd	s1,152(sp)
	mv	a0,sp
	.cfi_offset 9, -24
	lui	s1,%hi(__stack_chk_guard)
	ld	a5, %lo(__stack_chk_guard)(s1)
	sd	a5, 136(sp)
	li	a5, 0
	sd	ra,168(sp)
	.cfi_offset 1, -8
	call	fp_mul_4x4_u512_ct
	ld	a5,0(sp)
	addi	a1,sp,64
	mv	a0,s0
	sd	zero,128(sp)
	sd	a5,64(sp)
	ld	a5,8(sp)
	sd	a5,72(sp)
	ld	a5,16(sp)
	sd	a5,80(sp)
	ld	a5,24(sp)
	sd	a5,88(sp)
	ld	a5,32(sp)
	sd	a5,96(sp)
	ld	a5,40(sp)
	sd	a5,104(sp)
	ld	a5,48(sp)
	sd	a5,112(sp)
	ld	a5,56(sp)
	sd	a5,120(sp)
	call	fp_mont_reduce
	ld	a4, 136(sp)
	ld	a5, %lo(__stack_chk_guard)(s1)
	xor	a5, a4, a5
	li	a4, 0
	bne	a5,zero,.L212
	ld	ra,168(sp)
	.cfi_remember_state
	.cfi_restore 1
	ld	s0,160(sp)
	.cfi_restore 8
	ld	s1,152(sp)
	.cfi_restore 9
	addi	sp,sp,176
	.cfi_def_cfa_offset 0
	jr	ra
.L212:
	.cfi_restore_state
	call	__stack_chk_fail
	.cfi_endproc
.LFE26:
	.size	fp_mont_mul_ct, .-fp_mont_mul_ct
	.align	1
	.globl	fp_mont_mul
	.type	fp_mont_mul, @function
fp_mont_mul:
.LFB27:
	.cfi_startproc
	addi	sp,sp,-176
	.cfi_def_cfa_offset 176
	sd	s0,160(sp)
	.cfi_offset 8, -16
	mv	s0,a0
	sd	s1,152(sp)
	mv	a0,sp
	.cfi_offset 9, -24
	lui	s1,%hi(__stack_chk_guard)
	ld	a5, %lo(__stack_chk_guard)(s1)
	sd	a5, 136(sp)
	li	a5, 0
	sd	ra,168(sp)
	.cfi_offset 1, -8
	call	fp_mul_4x4_u512_ct
	ld	a5,0(sp)
	addi	a1,sp,64
	mv	a0,s0
	sd	zero,128(sp)
	sd	a5,64(sp)
	ld	a5,8(sp)
	sd	a5,72(sp)
	ld	a5,16(sp)
	sd	a5,80(sp)
	ld	a5,24(sp)
	sd	a5,88(sp)
	ld	a5,32(sp)
	sd	a5,96(sp)
	ld	a5,40(sp)
	sd	a5,104(sp)
	ld	a5,48(sp)
	sd	a5,112(sp)
	ld	a5,56(sp)
	sd	a5,120(sp)
	call	fp_mont_reduce
	ld	a4, 136(sp)
	ld	a5, %lo(__stack_chk_guard)(s1)
	xor	a5, a4, a5
	li	a4, 0
	bne	a5,zero,.L216
	ld	ra,168(sp)
	.cfi_remember_state
	.cfi_restore 1
	ld	s0,160(sp)
	.cfi_restore 8
	ld	s1,152(sp)
	.cfi_restore 9
	addi	sp,sp,176
	.cfi_def_cfa_offset 0
	jr	ra
.L216:
	.cfi_restore_state
	call	__stack_chk_fail
	.cfi_endproc
.LFE27:
	.size	fp_mont_mul, .-fp_mont_mul
	.align	1
	.globl	fp_mont_sqr
	.type	fp_mont_sqr, @function
fp_mont_sqr:
.LFB28:
	.cfi_startproc
	addi	sp,sp,-176
	.cfi_def_cfa_offset 176
	sd	s0,160(sp)
	.cfi_offset 8, -16
	mv	s0,a0
	sd	s1,152(sp)
	mv	a0,sp
	.cfi_offset 9, -24
	lui	s1,%hi(__stack_chk_guard)
	ld	a5, %lo(__stack_chk_guard)(s1)
	sd	a5, 136(sp)
	li	a5, 0
	sd	ra,168(sp)
	.cfi_offset 1, -8
	call	comba_sqr
	ld	a5,0(sp)
	addi	a1,sp,64
	mv	a0,s0
	sd	zero,128(sp)
	sd	a5,64(sp)
	ld	a5,8(sp)
	sd	a5,72(sp)
	ld	a5,16(sp)
	sd	a5,80(sp)
	ld	a5,24(sp)
	sd	a5,88(sp)
	ld	a5,32(sp)
	sd	a5,96(sp)
	ld	a5,40(sp)
	sd	a5,104(sp)
	ld	a5,48(sp)
	sd	a5,112(sp)
	ld	a5,56(sp)
	sd	a5,120(sp)
	call	fp_mont_reduce
	ld	a4, 136(sp)
	ld	a5, %lo(__stack_chk_guard)(s1)
	xor	a5, a4, a5
	li	a4, 0
	bne	a5,zero,.L220
	ld	ra,168(sp)
	.cfi_remember_state
	.cfi_restore 1
	ld	s0,160(sp)
	.cfi_restore 8
	ld	s1,152(sp)
	.cfi_restore 9
	addi	sp,sp,176
	.cfi_def_cfa_offset 0
	jr	ra
.L220:
	.cfi_restore_state
	call	__stack_chk_fail
	.cfi_endproc
.LFE28:
	.size	fp_mont_sqr, .-fp_mont_sqr
	.align	1
	.globl	fp_to_mont
	.type	fp_to_mont, @function
fp_to_mont:
.LFB29:
	.cfi_startproc
	addi	sp,sp,-176
	.cfi_def_cfa_offset 176
	lui	a2,%hi(.LANCHOR0+32)
	sd	s0,160(sp)
	addi	a2,a2,%lo(.LANCHOR0+32)
	sd	s1,152(sp)
	.cfi_offset 8, -16
	.cfi_offset 9, -24
	mv	s0,a0
	lui	s1,%hi(__stack_chk_guard)
	mv	a0,sp
	ld	a5, %lo(__stack_chk_guard)(s1)
	sd	a5, 136(sp)
	li	a5, 0
	sd	ra,168(sp)
	.cfi_offset 1, -8
	call	fp_mul_4x4_u512_ct
	sd	zero,128(sp)
	ld	a5,0(sp)
	addi	a1,sp,64
	mv	a0,s0
	sd	a5,64(sp)
	ld	a5,8(sp)
	sd	a5,72(sp)
	ld	a5,16(sp)
	sd	a5,80(sp)
	ld	a5,24(sp)
	sd	a5,88(sp)
	ld	a5,32(sp)
	sd	a5,96(sp)
	ld	a5,40(sp)
	sd	a5,104(sp)
	ld	a5,48(sp)
	sd	a5,112(sp)
	ld	a5,56(sp)
	sd	a5,120(sp)
	call	fp_mont_reduce
	ld	a4, 136(sp)
	ld	a5, %lo(__stack_chk_guard)(s1)
	xor	a5, a4, a5
	li	a4, 0
	bne	a5,zero,.L224
	ld	ra,168(sp)
	.cfi_remember_state
	.cfi_restore 1
	ld	s0,160(sp)
	.cfi_restore 8
	ld	s1,152(sp)
	.cfi_restore 9
	addi	sp,sp,176
	.cfi_def_cfa_offset 0
	jr	ra
.L224:
	.cfi_restore_state
	call	__stack_chk_fail
	.cfi_endproc
.LFE29:
	.size	fp_to_mont, .-fp_to_mont
	.align	1
	.globl	fp_from_mont
	.type	fp_from_mont, @function
fp_from_mont:
.LFB30:
	.cfi_startproc
	addi	sp,sp,-176
	.cfi_def_cfa_offset 176
	lui	a2,%hi(.LANCHOR0+64)
	sd	s0,160(sp)
	addi	a2,a2,%lo(.LANCHOR0+64)
	sd	s1,152(sp)
	.cfi_offset 8, -16
	.cfi_offset 9, -24
	mv	s0,a0
	lui	s1,%hi(__stack_chk_guard)
	mv	a0,sp
	ld	a5, %lo(__stack_chk_guard)(s1)
	sd	a5, 136(sp)
	li	a5, 0
	sd	ra,168(sp)
	.cfi_offset 1, -8
	call	fp_mul_4x4_u512_ct
	sd	zero,128(sp)
	ld	a5,0(sp)
	addi	a1,sp,64
	mv	a0,s0
	sd	a5,64(sp)
	ld	a5,8(sp)
	sd	a5,72(sp)
	ld	a5,16(sp)
	sd	a5,80(sp)
	ld	a5,24(sp)
	sd	a5,88(sp)
	ld	a5,32(sp)
	sd	a5,96(sp)
	ld	a5,40(sp)
	sd	a5,104(sp)
	ld	a5,48(sp)
	sd	a5,112(sp)
	ld	a5,56(sp)
	sd	a5,120(sp)
	call	fp_mont_reduce
	ld	a4, 136(sp)
	ld	a5, %lo(__stack_chk_guard)(s1)
	xor	a5, a4, a5
	li	a4, 0
	bne	a5,zero,.L228
	ld	ra,168(sp)
	.cfi_remember_state
	.cfi_restore 1
	ld	s0,160(sp)
	.cfi_restore 8
	ld	s1,152(sp)
	.cfi_restore 9
	addi	sp,sp,176
	.cfi_def_cfa_offset 0
	jr	ra
.L228:
	.cfi_restore_state
	call	__stack_chk_fail
	.cfi_endproc
.LFE30:
	.size	fp_from_mont, .-fp_from_mont
	.align	1
	.globl	fp_mont_inv
	.type	fp_mont_inv, @function
fp_mont_inv:
.LFB31:
	.cfi_startproc
	addi	sp,sp,-576
	.cfi_def_cfa_offset 576
	sd	s2,544(sp)
	.cfi_offset 18, -32
	lui	s2,%hi(__stack_chk_guard)
	sd	s3,536(sp)
	ld	a5, %lo(__stack_chk_guard)(s2)
	sd	a5, 520(sp)
	li	a5, 0
	.cfi_offset 19, -40
	mv	s3,a0
	addi	a0,sp,384
	sd	ra,568(sp)
	sd	s0,560(sp)
	.cfi_offset 1, -8
	.cfi_offset 8, -16
	li	s0,3
	sd	s1,552(sp)
	.cfi_offset 9, -24
	mv	s1,a1
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	mv	a0,sp
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	mv	a2,s1
	mv	a1,sp
	addi	a0,sp,384
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	mv	a0,sp
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	mv	a1,sp
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,32
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,32
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,32
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	mv	a2,sp
	addi	a1,sp,32
	addi	a0,sp,384
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,32
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,32
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,64
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L230:
	addi	a1,sp,64
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,64
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L230
	addi	a2,sp,32
	addi	a1,sp,64
	addi	a0,sp,384
	li	s0,7
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,64
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,64
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,96
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L231:
	addi	a1,sp,96
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,96
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L231
	addi	a2,sp,64
	addi	a1,sp,96
	addi	a0,sp,384
	li	s0,15
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,96
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,96
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,128
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L232:
	addi	a1,sp,128
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,128
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L232
	addi	a2,sp,96
	addi	a1,sp,128
	addi	a0,sp,384
	li	s0,5
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,128
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,32
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,160
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,160
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,160
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	mv	a2,sp
	addi	a1,sp,160
	addi	a0,sp,384
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,160
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,64
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,192
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L233:
	addi	a1,sp,192
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,192
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L233
	addi	a2,sp,160
	addi	a1,sp,192
	addi	a0,sp,384
	li	s0,13
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,192
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,96
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,224
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L234:
	addi	a1,sp,224
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,224
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L234
	addi	a2,sp,192
	addi	a1,sp,224
	addi	a0,sp,384
	li	s0,31
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,224
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,224
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,256
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	mv	a2,s1
	addi	a1,sp,256
	addi	a0,sp,384
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,256
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,128
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,288
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L235:
	addi	a1,sp,288
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,288
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L235
	addi	a2,sp,128
	addi	a1,sp,288
	addi	a0,sp,384
	li	s0,32
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,288
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,256
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,352
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L236:
	addi	a1,sp,352
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,352
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L236
	addi	a2,sp,352
	addi	a1,sp,128
	addi	a0,sp,384
	li	s0,63
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L237:
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L237
	addi	a2,sp,288
	addi	a1,sp,320
	addi	a0,sp,384
	li	s0,31
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L238:
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L238
	addi	a2,sp,128
	addi	a1,sp,320
	addi	a0,sp,384
	li	s0,31
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L239:
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L239
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	li	s0,31
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L240:
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L240
	addi	a2,sp,128
	addi	a1,sp,320
	addi	a0,sp,384
	li	s0,29
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
.L241:
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	addiw	s0,s0,-1
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	zero,512(sp)
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	bne	s0,zero,.L241
	addi	a2,sp,224
	addi	a1,sp,320
	addi	a0,sp,384
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	addi	a1,sp,320
	addi	a0,sp,384
	call	comba_sqr
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	addi	a0,sp,320
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	mv	a2,s1
	addi	a1,sp,320
	addi	a0,sp,384
	call	fp_mul_4x4_u512_ct
	sd	zero,512(sp)
	ld	a5,384(sp)
	addi	a1,sp,448
	mv	a0,s3
	sd	a5,448(sp)
	ld	a5,392(sp)
	sd	a5,456(sp)
	ld	a5,400(sp)
	sd	a5,464(sp)
	ld	a5,408(sp)
	sd	a5,472(sp)
	ld	a5,416(sp)
	sd	a5,480(sp)
	ld	a5,424(sp)
	sd	a5,488(sp)
	ld	a5,432(sp)
	sd	a5,496(sp)
	ld	a5,440(sp)
	sd	a5,504(sp)
	call	fp_mont_reduce
	ld	a4, 520(sp)
	ld	a5, %lo(__stack_chk_guard)(s2)
	xor	a5, a4, a5
	li	a4, 0
	bne	a5,zero,.L256
	ld	ra,568(sp)
	.cfi_remember_state
	.cfi_restore 1
	ld	s0,560(sp)
	.cfi_restore 8
	ld	s1,552(sp)
	.cfi_restore 9
	ld	s2,544(sp)
	.cfi_restore 18
	ld	s3,536(sp)
	.cfi_restore 19
	addi	sp,sp,576
	.cfi_def_cfa_offset 0
	jr	ra
.L256:
	.cfi_restore_state
	call	__stack_chk_fail
	.cfi_endproc
.LFE31:
	.size	fp_mont_inv, .-fp_mont_inv
	.globl	FP_MONT_R2
	.globl	FP_MONT_ONE
	.globl	FP_ONE
	.globl	FP_ZERO
	.globl	FP_P
	.section	.rodata
	.align	3
	.set	.LANCHOR0,. + 0
	.type	FP_MONT_ONE, @object
	.size	FP_MONT_ONE, 32
FP_MONT_ONE:
	.dword	1
	.dword	4294967295
	.dword	0
	.dword	4294967296
	.type	FP_MONT_R2, @object
	.size	FP_MONT_R2, 32
FP_MONT_R2:
	.dword	8589934595
	.dword	12884901887
	.dword	4294967297
	.dword	17179869186
	.type	FP_ONE, @object
	.size	FP_ONE, 32
FP_ONE:
	.dword	1
	.dword	0
	.dword	0
	.dword	0
	.type	FP_ZERO, @object
	.size	FP_ZERO, 32
FP_ZERO:
	.zero	32
	.type	FP_P, @object
	.size	FP_P, 32
FP_P:
	.dword	-1
	.dword	-4294967296
	.dword	-1
	.dword	-4294967297
	.ident	"GCC: (Bianbu 13.2.0-23ubuntu4bb3) 13.2.0"
	.section	.note.GNU-stack,"",@progbits
