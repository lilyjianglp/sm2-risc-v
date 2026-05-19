# =================================================================
# SM2 RISC-V 64 Optimization Project Makefile
# Muse Pi Pro / SpacemiT K1(X60) board-only tuned version
#
# fn 层统一固定为：fn_Montgomery.c + asm/fn_mont.S
# 这样可以保证 fn_mont_mul_ct 等符号由 asm/fn_mont.S 提供
# =================================================================

CC      ?= gcc
CXX     ?= g++
AR      ?= ar
SIM     :=
SYSROOT :=

USE_STATIC ?= 0
USE_LTO    ?= 0
USE_PGO    ?= 0

# -----------------------------------------------------------------
# Board-specific CPU tuning
# 固定为 SpacemiT X60，不再保留 native / generic 分支
# -----------------------------------------------------------------
CPU_TUNE   := spacemit-x60
CPU_CFLAGS := -mcpu=spacemit-x60

# -----------------------------------------------------------------
# Optimization / common flags
# -----------------------------------------------------------------
BASE_CFLAGS = -Wall -Wextra -Wshadow -Wno-unused-parameter -std=gnu11
OPT_CFLAGS  = -O3 -fomit-frame-pointer -fno-plt -falign-functions=32

# 对这种大整数/密码热点代码，通常不建议开严格别名假设，减少潜在误判
SAFE_CFLAGS = -fno-strict-aliasing

# 汇编文件走预处理，方便 #ifdef / 宏展开
ASM_CPPFLAGS = -x assembler-with-cpp

ifeq ($(USE_LTO),1)
	LTO_FLAGS = -flto
else
	LTO_FLAGS =
endif

# PGO 两阶段：
# 1) make clean && make bench-all USE_PGO=gen
# 2) 跑 ./test_xxx 生成 profile
# 3) make clean && make bench-all USE_PGO=use
ifeq ($(USE_PGO),gen)
	PGO_FLAGS = -fprofile-generate
else ifeq ($(USE_PGO),use)
	PGO_FLAGS = -fprofile-use -fprofile-correction
else
	PGO_FLAGS =
endif

# -----------------------------------------------------------------
# 宏定义统一走 CPPFLAGS
# 这样 C 和 .S 都能更清楚地继承
# -----------------------------------------------------------------
CPPFLAGS =

# C 编译参数
CFLAGS = \
	$(BASE_CFLAGS) \
	$(OPT_CFLAGS) \
	$(SAFE_CFLAGS) \
	$(CPU_CFLAGS) \
	$(LTO_FLAGS) \
	$(PGO_FLAGS)

CXXFLAGS = \
	-Wall -Wextra -Wshadow -std=c++17 \
	$(OPT_CFLAGS) \
	$(SAFE_CFLAGS) \
	$(CPU_CFLAGS) \
	$(LTO_FLAGS) \
	$(PGO_FLAGS) \
	-I.

PROFILE_CFLAGS = -g -fno-omit-frame-pointer -fno-inline

# ASM 编译参数
# 注意：不给 .S 传 -std=gnu11 / -Wshadow 这种纯 C 警告参数
# 但要传优化、CPU、LTO/PGO、assembler-with-cpp
ASFLAGS = \
	$(OPT_CFLAGS) \
	$(SAFE_CFLAGS) \
	$(CPU_CFLAGS) \
	$(LTO_FLAGS) \
	$(PGO_FLAGS) \
	$(ASM_CPPFLAGS)

ifeq ($(USE_STATIC),1)
	LINK_MODE = -static
else
	LINK_MODE =
endif

LDFLAGS = $(LINK_MODE) $(LTO_FLAGS) $(PGO_FLAGS)

QEMU_RUN =

# -------------------------------------------------
# 测试入口
# -------------------------------------------------
TEST_SRC = test.c

# -------------------------------------------------
# fn 层：所有 target 全部固定
# fn_Montgomery.c + asm/fn_mont.S
# -------------------------------------------------
FN_FIXED_C      = fn_Montgomery.c
FN_FIXED_ASM    = asm/fn_mont.S
FN_FIXED_DEFS   = -DHAVE_FN_MONT_CORE -DUSE_ASM_FN_MONT -DNO_TEST_FN_INTERNALS

# -------------------------------------------------
# Prime / 普通模 p 路线
# -------------------------------------------------
FP_PRIME_C     = fp.c
CURVE_PRIME_C  = sm2_curve.c
SCALAR_PRIME_C = sm2_scalar.c
KEX_PRIME_C    = sm2_kex.c

PRIME_SRCS = \
	$(FP_PRIME_C) \
	$(CURVE_PRIME_C) \
	$(SCALAR_PRIME_C) \
	$(KEX_PRIME_C) \
	$(FN_FIXED_C) \
	$(TEST_SRC)

# -------------------------------------------------
# Montgomery 模 p 路线：无 Pornin
# 这里不要加 -DUSE_PORNIN_INV，也不要链接 fp_pornin_inv.c
# 使用 fp_montgomery.c 里 #ifndef USE_PORNIN_INV 包住的原始 fp_mont_inv()
# -------------------------------------------------
FP_MONT_C      = fp_montgomery.c
CURVE_MONT_C   = sm2_curve_mont.c
SCALAR_MONT_C  = sm2_scalar_mont.c
KEX_MONT_C     = sm2_kex_mont.c

MONT_SRCS = \
	$(FP_MONT_C) \
	$(CURVE_MONT_C) \
	$(SCALAR_MONT_C) \
	$(KEX_MONT_C) \
	$(FN_FIXED_C) \
	$(TEST_SRC)

# -------------------------------------------------
# Montgomery + Pornin 路线
# 这里必须加 -DUSE_PORNIN_INV，并链接 fp_pornin_inv.c
# 使用 fp_pornin_inv.c 里的 Pornin fp_mont_inv()
# -------------------------------------------------
FP_PORNIN_C = fp_pornin_inv.c

# Pornin 31-step inner loop ASM
PORNIN_INNER31_ASM = asm/pornin_inner31.S
PORNIN_INNER31_DEFS = -DUSE_PORNIN_INNER31_ASM

# Inner31 in the full inversion asm can be emitted either fully unrolled
# or as a smaller 4-step partially unrolled loop. Keep full unroll as the
# default; enable this only for front-end/I-cache A/B testing:
#   make perf-pornin-asm USE_FULL_INV_ASM=1 USE_INNER31_HALF_UNROLL=1
USE_INNER31_HALF_UNROLL ?= 0
ifeq ($(USE_INNER31_HALF_UNROLL),1)
PORNIN_INNER31_UNROLL_DEFS = -DUSE_PORNIN_INNER31_HALF_UNROLL
else
PORNIN_INNER31_UNROLL_DEFS =
endif

# Pornin small signed lincomb ASM, used by update_uv.
PORNIN_LINCOMB_ASM = asm/pornin_lincomb.S
PORNIN_LINCOMB_DEFS = -DUSE_PORNIN_LINCOMB_ASM
PORNIN_UPDATE_UV_DEFS = -DUSE_PORNIN_UPDATE_UV_ASM
PORNIN_UPDATE_AB_ASM = asm/pornin_update_ab.S
PORNIN_UPDATE_AB_DEFS = -DUSE_PORNIN_UPDATE_AB_ASM
PORNIN_FULL_INV_ASM = asm/pornin_full_inv.S
PORNIN_FULL_INV_DEFS = -DUSE_PORNIN_FULL_INV_ASM

# Experimental update_ab ASM. It is correct, but slower than the C version
# on X60 so keep it disabled by default.
# Enable only for A/B testing:
#   make perf-pornin-asm USE_AB_ASM=1
USE_AB_ASM ?= 0
ifeq ($(USE_AB_ASM),1)
PORNIN_UPDATE_AB_SRCS = $(PORNIN_UPDATE_AB_ASM)
PORNIN_UPDATE_AB_ACTIVE_DEFS = $(PORNIN_UPDATE_AB_DEFS)
else
PORNIN_UPDATE_AB_SRCS =
PORNIN_UPDATE_AB_ACTIVE_DEFS =
endif

# Experimental full fp_mont_inv() ASM driver.
# NOTE: the previous helper-calling driver is intentionally not part of the
# default path. The real target should inline extract/inner31/update_ab/update_uv
# instead of calling update_ab_asm/update_uv_asm.
# Enable only when working on that fused file:
#   make run-pornin-inv-asm USE_FULL_INV_ASM=1
USE_FULL_INV_ASM ?= 0
ifeq ($(USE_FULL_INV_ASM),1)
PORNIN_FULL_INV_SRCS = $(PORNIN_FULL_INV_ASM)
PORNIN_FULL_INV_ACTIVE_DEFS = $(PORNIN_FULL_INV_DEFS)
else
PORNIN_FULL_INV_SRCS =
PORNIN_FULL_INV_ACTIVE_DEFS =
endif

PORNIN_EXTRA_ASM_SRCS = $(sort $(PORNIN_UPDATE_AB_SRCS) $(PORNIN_FULL_INV_SRCS))

MONT_PORNIN_SRCS = \
	$(FP_MONT_C) \
	$(FP_PORNIN_C) \
	$(CURVE_MONT_C) \
	$(SCALAR_MONT_C) \
	$(KEX_MONT_C) \
	$(FN_FIXED_C) \
	$(TEST_SRC)

# -------------------------------------------------
# Prime fp ASM
# -------------------------------------------------
ASM_FP_PRIME_SRCS = \
	asm/fp_add.S \
	asm/fp_sub.S \
	asm/fp_neg.S \
	asm/fp_mul.S \
	asm/fp_sqr.S \
	asm/fp_reduce.S

ASM_FP_PRIME_DEFS = \
	-DUSE_ASM_FP \
	-DUSE_ASM_REDUCE

# -------------------------------------------------
# Montgomery fp ASM
# -------------------------------------------------
ASM_FP_MONT_SRCS = \
	asm/fp_mont_add.S \
	asm/fp_mont_sub.S \
	asm/fp_mont_neg.S \
	asm/fp_mont_mul.S \
	asm/fp_mont_sqr.S \
	asm/fp_mont_to_from.S

ASM_FP_MONT_DEFS = \
	-DUSE_RV64_ASM_ADD \
	-DUSE_RV64_ASM_SUB \
	-DUSE_RV64_ASM_NEG \
	-DUSE_RV64_ASM_MUL \
	-DUSE_RV64_ASM_SQR \
	-DUSE_RV64_ASM_TO_FROM

# -------------------------------------------------
# Targets
# -------------------------------------------------
TARGET_PRIME_C       = test_prime_c
TARGET_PRIME_FP_ASM  = test_prime_fp_asm
TARGET_MONT_C        = test_mont_c
TARGET_MONT_ASM      = test_mont_asm
TARGET_PORNIN_C      = test_pornin_c
TARGET_PORNIN_ASM    = test_pornin_asm
TARGET_PORNIN_INV_C   = test_pornin_inv_c
TARGET_PORNIN_INV_ASM = test_pornin_inv_asm
TARGET_PORNIN_INV_PROF_C   = test_pornin_inv_prof_c
TARGET_PORNIN_INV_PROF_ASM = test_pornin_inv_prof_asm
TARGET_INV_PERF = test_inv_perf
TARGET_INV_PERF_MONT_C = test_inv_perf_mont_c
TARGET_INV_PERF_PORNIN_C = test_inv_perf_pornin_c
TARGET_SECURE_MESSAGE_DEMO = secure_message_demo

SECURE_BUILD_DIR = build/secure_message
SECURE_C_OBJS = \
	$(SECURE_BUILD_DIR)/fp_montgomery.o \
	$(SECURE_BUILD_DIR)/fp_pornin_inv.o \
	$(SECURE_BUILD_DIR)/fn_Montgomery.o \
	$(SECURE_BUILD_DIR)/sm2_curve_mont.o \
	$(SECURE_BUILD_DIR)/sm2_scalar_mont.o \
	$(SECURE_BUILD_DIR)/sm2_kex_mont.o

SECURE_ASM_OBJS = \
	$(SECURE_BUILD_DIR)/fp_mont_add.o \
	$(SECURE_BUILD_DIR)/fp_mont_sub.o \
	$(SECURE_BUILD_DIR)/fp_mont_neg.o \
	$(SECURE_BUILD_DIR)/fp_mont_mul.o \
	$(SECURE_BUILD_DIR)/fp_mont_sqr.o \
	$(SECURE_BUILD_DIR)/fp_mont_to_from.o \
	$(SECURE_BUILD_DIR)/pornin_full_inv.o \
	$(SECURE_BUILD_DIR)/fn_mont.o

SECURE_CPP_OBJS = \
	$(SECURE_BUILD_DIR)/sm4_channel.o \
	$(SECURE_BUILD_DIR)/secure_message_demo.o

all: $(TARGET_PRIME_C) $(TARGET_PRIME_FP_ASM) \
     $(TARGET_MONT_C) $(TARGET_MONT_ASM) \
     $(TARGET_PORNIN_C) $(TARGET_PORNIN_ASM)

# -------------------------------------------------
# Prime 纯 C（fn 固定为 C+asm）
# -------------------------------------------------
$(TARGET_PRIME_C): $(PRIME_SRCS) $(FN_FIXED_ASM)
	@echo ">>> Building Prime / Pure C (fn fixed to C+asm) ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      $(FN_FIXED_DEFS) \
	      $(PRIME_SRCS) $(FN_FIXED_ASM) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# Prime fp ASM（fn 固定为 C+asm）
# -------------------------------------------------
$(TARGET_PRIME_FP_ASM): $(PRIME_SRCS) $(ASM_FP_PRIME_SRCS) $(FN_FIXED_ASM)
	@echo ">>> Building Prime / fp ASM (fn fixed to C+asm) ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      $(FN_FIXED_DEFS) \
	      $(ASM_FP_PRIME_DEFS) \
	      $(PRIME_SRCS) \
	      $(ASFLAGS) $(ASM_FP_PRIME_SRCS) $(FN_FIXED_ASM) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# Montgomery 纯 C：无 Pornin
# -------------------------------------------------
$(TARGET_MONT_C): $(MONT_SRCS) $(FN_FIXED_ASM)
	@echo ">>> Building Montgomery / Pure C, no Pornin (fn fixed to C+asm) ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT $(FN_FIXED_DEFS) \
	      $(MONT_SRCS) $(FN_FIXED_ASM) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# Montgomery ASM：无 Pornin
# -------------------------------------------------
$(TARGET_MONT_ASM): $(MONT_SRCS) $(ASM_FP_MONT_SRCS) $(FN_FIXED_ASM)
	@echo ">>> Building Montgomery / fp ASM, no Pornin (fn fixed to C+asm) ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT $(FN_FIXED_DEFS) \
	      $(ASM_FP_MONT_DEFS) \
	      $(MONT_SRCS) \
	      $(ASFLAGS) $(ASM_FP_MONT_SRCS) $(FN_FIXED_ASM) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# Pornin 纯 C：Montgomery + Pornin inv
# -------------------------------------------------
$(TARGET_PORNIN_C): $(MONT_PORNIN_SRCS) $(FN_FIXED_ASM)
	@echo ">>> Building Montgomery / Pure C, Pornin inv (fn fixed to C+asm) ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV $(FN_FIXED_DEFS) \
	      $(MONT_PORNIN_SRCS) $(FN_FIXED_ASM) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# Pornin ASM：Montgomery + Pornin inv + fp ASM + inner31 ASM
# -------------------------------------------------
$(TARGET_PORNIN_ASM): $(MONT_PORNIN_SRCS) $(ASM_FP_MONT_SRCS) $(FN_FIXED_ASM) $(PORNIN_INNER31_ASM) $(PORNIN_LINCOMB_ASM) $(PORNIN_EXTRA_ASM_SRCS)
	@echo ">>> Building Montgomery / fp ASM, Pornin inv + inner31 ASM (fn fixed to C+asm) ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV $(PORNIN_INNER31_DEFS) $(PORNIN_LINCOMB_DEFS) $(PORNIN_UPDATE_UV_DEFS) $(PORNIN_UPDATE_AB_ACTIVE_DEFS) $(PORNIN_FULL_INV_ACTIVE_DEFS) $(FN_FIXED_DEFS) \
	      $(PORNIN_INNER31_UNROLL_DEFS) \
	      $(ASM_FP_MONT_DEFS) \
	      $(MONT_PORNIN_SRCS) \
	      $(ASFLAGS) $(ASM_FP_MONT_SRCS) $(PORNIN_INNER31_ASM) $(PORNIN_LINCOMB_ASM) $(PORNIN_EXTRA_ASM_SRCS) $(FN_FIXED_ASM) \
	      $(LDFLAGS) -o $@

$(TARGET_PORNIN_INV_C): test_pornin_inv.c $(FP_MONT_C) $(FP_PORNIN_C)
	@echo ">>> Building Pornin inversion correctness / C inner loop ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV \
	      test_pornin_inv.c $(FP_MONT_C) $(FP_PORNIN_C) \
	      $(LDFLAGS) -o $@

$(TARGET_PORNIN_INV_ASM): test_pornin_inv.c $(FP_MONT_C) $(FP_PORNIN_C) $(ASM_FP_MONT_SRCS) $(PORNIN_INNER31_ASM) $(PORNIN_LINCOMB_ASM) $(PORNIN_EXTRA_ASM_SRCS)
	@echo ">>> Building Pornin inversion correctness / ASM inner31 ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV $(PORNIN_INNER31_DEFS) $(PORNIN_LINCOMB_DEFS) $(PORNIN_UPDATE_UV_DEFS) $(PORNIN_UPDATE_AB_ACTIVE_DEFS) $(PORNIN_FULL_INV_ACTIVE_DEFS) \
	      $(PORNIN_INNER31_UNROLL_DEFS) \
	      $(ASM_FP_MONT_DEFS) \
	      test_pornin_inv.c $(FP_MONT_C) $(FP_PORNIN_C) \
	      $(ASFLAGS) $(ASM_FP_MONT_SRCS) $(PORNIN_INNER31_ASM) $(PORNIN_LINCOMB_ASM) $(PORNIN_EXTRA_ASM_SRCS) \
	      $(LDFLAGS) -o $@

$(TARGET_PORNIN_INV_PROF_C): test_pornin_inv.c $(FP_MONT_C) $(FP_PORNIN_C)
	@echo ">>> Building Pornin inversion profile / C inner loop ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) $(PROFILE_CFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV -DPORNIN_PROFILE -DPORNIN_PERF_NOINLINE \
	      test_pornin_inv.c $(FP_MONT_C) $(FP_PORNIN_C) \
	      $(LDFLAGS) -o $@

$(TARGET_PORNIN_INV_PROF_ASM): test_pornin_inv.c $(FP_MONT_C) $(FP_PORNIN_C) $(ASM_FP_MONT_SRCS) $(PORNIN_INNER31_ASM) $(PORNIN_LINCOMB_ASM) $(PORNIN_EXTRA_ASM_SRCS)
	@echo ">>> Building Pornin inversion profile / ASM inner31 ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) $(PROFILE_CFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV -DPORNIN_PROFILE -DPORNIN_PERF_NOINLINE $(PORNIN_INNER31_DEFS) $(PORNIN_LINCOMB_DEFS) $(PORNIN_UPDATE_UV_DEFS) $(PORNIN_UPDATE_AB_ACTIVE_DEFS) $(PORNIN_FULL_INV_ACTIVE_DEFS) \
	      $(PORNIN_INNER31_UNROLL_DEFS) \
	      $(ASM_FP_MONT_DEFS) \
	      test_pornin_inv.c $(FP_MONT_C) $(FP_PORNIN_C) \
	      $(ASFLAGS) $(ASM_FP_MONT_SRCS) $(PORNIN_INNER31_ASM) $(PORNIN_LINCOMB_ASM) $(PORNIN_EXTRA_ASM_SRCS) \
	      $(LDFLAGS) -o $@

$(TARGET_INV_PERF): test_inv_perf.c $(FP_MONT_C) $(FP_PORNIN_C) $(ASM_FP_MONT_SRCS) $(PORNIN_INNER31_ASM) $(PORNIN_LINCOMB_ASM) $(PORNIN_EXTRA_ASM_SRCS)
	@echo ">>> Building focused fp_mont_inv benchmark ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV $(PORNIN_INNER31_DEFS) $(PORNIN_LINCOMB_DEFS) $(PORNIN_UPDATE_UV_DEFS) $(PORNIN_UPDATE_AB_ACTIVE_DEFS) $(PORNIN_FULL_INV_ACTIVE_DEFS) \
	      $(PORNIN_INNER31_UNROLL_DEFS) \
	      $(ASM_FP_MONT_DEFS) \
	      test_inv_perf.c $(FP_MONT_C) $(FP_PORNIN_C) \
	      $(ASFLAGS) $(ASM_FP_MONT_SRCS) $(PORNIN_INNER31_ASM) $(PORNIN_LINCOMB_ASM) $(PORNIN_EXTRA_ASM_SRCS) \
	      $(LDFLAGS) -o $@

$(TARGET_INV_PERF_MONT_C): test_inv_perf.c $(FP_MONT_C)
	@echo ">>> Building focused fp_mont_inv benchmark / no Pornin C ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT \
	      test_inv_perf.c $(FP_MONT_C) \
	      $(LDFLAGS) -o $@

$(TARGET_INV_PERF_PORNIN_C): test_inv_perf.c $(FP_MONT_C) $(FP_PORNIN_C)
	@echo ">>> Building focused fp_mont_inv benchmark / Pornin C ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV \
	      test_inv_perf.c $(FP_MONT_C) $(FP_PORNIN_C) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# test_perf：默认保留原来的 Montgomery / fp ASM 无 Pornin 版本
# -------------------------------------------------
TARGET_PERF_MONT_ASM = test_perf
MONT_PERF_SRCS = $(filter-out $(TEST_SRC),$(MONT_SRCS))

$(TARGET_PERF_MONT_ASM): test_perf.c $(MONT_PERF_SRCS) $(ASM_FP_MONT_SRCS) $(FN_FIXED_ASM)
	@echo ">>> Building test_perf (Montgomery / fp ASM, no Pornin, fn fixed to C+asm) ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT $(FN_FIXED_DEFS) \
	      $(ASM_FP_MONT_DEFS) \
	      test_perf.c \
	      $(MONT_PERF_SRCS) \
	      $(ASFLAGS) $(ASM_FP_MONT_SRCS) $(FN_FIXED_ASM) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# SM4 secure message demo
# -------------------------------------------------
$(SECURE_BUILD_DIR):
	mkdir -p $(SECURE_BUILD_DIR)

$(SECURE_BUILD_DIR)/%.o: %.c | $(SECURE_BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV $(PORNIN_FULL_INV_DEFS) \
	      $(FN_FIXED_DEFS) $(ASM_FP_MONT_DEFS) \
	      -c $< -o $@

$(SECURE_BUILD_DIR)/%.o: asm/%.S | $(SECURE_BUILD_DIR)
	$(CC) $(CPPFLAGS) $(ASFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV $(PORNIN_FULL_INV_DEFS) \
	      $(FN_FIXED_DEFS) $(ASM_FP_MONT_DEFS) \
	      -c $< -o $@

$(SECURE_BUILD_DIR)/sm4_channel.o: sm4_channel.cpp sm4_channel.h | $(SECURE_BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV $(PORNIN_FULL_INV_DEFS) \
	      $(FN_FIXED_DEFS) $(ASM_FP_MONT_DEFS) \
	      -c $< -o $@

$(SECURE_BUILD_DIR)/secure_message_demo.o: demo/secure_message_demo.cpp sm4_channel.h | $(SECURE_BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) \
	      -DUSE_FP_MONT -DUSE_PORNIN_INV $(PORNIN_FULL_INV_DEFS) \
	      $(FN_FIXED_DEFS) $(ASM_FP_MONT_DEFS) \
	      -c $< -o $@

$(TARGET_SECURE_MESSAGE_DEMO): $(SECURE_C_OBJS) $(SECURE_ASM_OBJS) $(SECURE_CPP_OBJS)
	@echo ">>> Building SM4 secure message demo (Montgomery / fp ASM / Pornin full inv ASM) ..."
	$(CXX) $(LDFLAGS) $^ -o $@

# -------------------------------------------------
# Run
# -------------------------------------------------
run-prime-c: $(TARGET_PRIME_C)
	./$(TARGET_PRIME_C)

run-prime-asm: $(TARGET_PRIME_FP_ASM)
	./$(TARGET_PRIME_FP_ASM)

run-mont-c: $(TARGET_MONT_C)
	./$(TARGET_MONT_C)

run-mont-asm: $(TARGET_MONT_ASM)
	./$(TARGET_MONT_ASM)

run-pornin-c: $(TARGET_PORNIN_C)
	./$(TARGET_PORNIN_C)

run-pornin-asm: $(TARGET_PORNIN_ASM)
	./$(TARGET_PORNIN_ASM)

run-pornin-inv-c: $(TARGET_PORNIN_INV_C)
	./$(TARGET_PORNIN_INV_C)

run-pornin-inv-asm: $(TARGET_PORNIN_INV_ASM)
	./$(TARGET_PORNIN_INV_ASM)

run-pornin-inv-prof-c: $(TARGET_PORNIN_INV_PROF_C)
	./$(TARGET_PORNIN_INV_PROF_C)

run-pornin-inv-prof-asm: $(TARGET_PORNIN_INV_PROF_ASM)
	./$(TARGET_PORNIN_INV_PROF_ASM)

run-test-perf: $(TARGET_PERF_MONT_ASM)
	./$(TARGET_PERF_MONT_ASM)

run-secure-message-demo: $(TARGET_SECURE_MESSAGE_DEMO)
	./$(TARGET_SECURE_MESSAGE_DEMO)

run-inv-perf: $(TARGET_INV_PERF)
	$(PERF) ./$(TARGET_INV_PERF)

run-inv-perf-mont-c: $(TARGET_INV_PERF_MONT_C)
	$(PERF) ./$(TARGET_INV_PERF_MONT_C)

run-inv-perf-pornin-c: $(TARGET_INV_PERF_PORNIN_C)
	$(PERF) ./$(TARGET_INV_PERF_PORNIN_C)

# -------------------------------------------------
# Bench
# -------------------------------------------------
bench-prime: $(TARGET_PRIME_C) $(TARGET_PRIME_FP_ASM)
	@echo "========================================================"
	@echo "CPU_TUNE=$(CPU_TUNE)  fn=C+asm fixed"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "ASFLAGS=$(ASFLAGS)"
	@echo "=== Prime / Pure C ==="
	./$(TARGET_PRIME_C)
	@echo ""
	@echo "=== Prime / fp ASM ==="
	./$(TARGET_PRIME_FP_ASM)
	@echo "========================================================"

bench-mont: $(TARGET_MONT_C) $(TARGET_MONT_ASM) $(TARGET_PORNIN_C) $(TARGET_PORNIN_ASM)
	@echo "========================================================"
	@echo "CPU_TUNE=$(CPU_TUNE)  fn=C+asm fixed"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "ASFLAGS=$(ASFLAGS)"
	@echo "=== Montgomery / Pure C, no Pornin ==="
	./$(TARGET_MONT_C)
	@echo ""
	@echo "=== Montgomery / Pure C, Pornin ==="
	./$(TARGET_PORNIN_C)
	@echo ""
	@echo "=== Montgomery / fp ASM, no Pornin ==="
	./$(TARGET_MONT_ASM)
	@echo ""
	@echo "=== Montgomery / fp ASM, Pornin ==="
	./$(TARGET_PORNIN_ASM)
	@echo "========================================================"

bench-all: $(TARGET_PRIME_C) $(TARGET_MONT_C) $(TARGET_PORNIN_C) \
           $(TARGET_PRIME_FP_ASM) $(TARGET_MONT_ASM) $(TARGET_PORNIN_ASM)
	@echo "========================================================"
	@echo "CPU_TUNE=$(CPU_TUNE)  fn=C+asm fixed"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "ASFLAGS=$(ASFLAGS)"
	@echo "=== Prime / Pure C ==="
	./$(TARGET_PRIME_C)
	@echo ""
	@echo "=== Montgomery / Pure C, no Pornin ==="
	./$(TARGET_MONT_C)
	@echo ""
	@echo "=== Montgomery / Pure C, Pornin ==="
	./$(TARGET_PORNIN_C)
	@echo ""
	@echo "=== Prime / fp ASM ==="
	./$(TARGET_PRIME_FP_ASM)
	@echo ""
	@echo "=== Montgomery / fp ASM, no Pornin ==="
	./$(TARGET_MONT_ASM)
	@echo ""
	@echo "=== Montgomery / fp ASM, Pornin ==="
	./$(TARGET_PORNIN_ASM)
	@echo "========================================================"

# -------------------------------------------------
# Perf
# -------------------------------------------------
PERF = perf stat -e cycles,instructions

perf-prime-c: $(TARGET_PRIME_C)
	$(PERF) ./$(TARGET_PRIME_C)

perf-mont-c: $(TARGET_MONT_C)
	$(PERF) ./$(TARGET_MONT_C)

perf-pornin-c: $(TARGET_PORNIN_C)
	$(PERF) ./$(TARGET_PORNIN_C)

perf-prime-asm: $(TARGET_PRIME_FP_ASM)
	$(PERF) ./$(TARGET_PRIME_FP_ASM)

perf-mont-asm: $(TARGET_MONT_ASM)
	$(PERF) ./$(TARGET_MONT_ASM)

perf-pornin-asm: $(TARGET_PORNIN_ASM)
	$(PERF) ./$(TARGET_PORNIN_ASM)

perf-inv: $(TARGET_INV_PERF)
	$(PERF) ./$(TARGET_INV_PERF)

perf-inv-mont-c: $(TARGET_INV_PERF_MONT_C)
	$(PERF) ./$(TARGET_INV_PERF_MONT_C)

perf-inv-pornin-c: $(TARGET_INV_PERF_PORNIN_C)
	$(PERF) ./$(TARGET_INV_PERF_PORNIN_C)

perf-all: $(TARGET_PRIME_C) $(TARGET_MONT_C) $(TARGET_PORNIN_C) \
          $(TARGET_PRIME_FP_ASM) $(TARGET_MONT_ASM) $(TARGET_PORNIN_ASM)
	@echo "=== Prime / Pure C ==="
	$(PERF) ./$(TARGET_PRIME_C)
	@echo ""
	@echo "=== Montgomery / Pure C, no Pornin ==="
	$(PERF) ./$(TARGET_MONT_C)
	@echo ""
	@echo "=== Montgomery / Pure C, Pornin ==="
	$(PERF) ./$(TARGET_PORNIN_C)
	@echo ""
	@echo "=== Prime / fp ASM ==="
	$(PERF) ./$(TARGET_PRIME_FP_ASM)
	@echo ""
	@echo "=== Montgomery / fp ASM, no Pornin ==="
	$(PERF) ./$(TARGET_MONT_ASM)
	@echo ""
	@echo "=== Montgomery / fp ASM, Pornin ==="
	$(PERF) ./$(TARGET_PORNIN_ASM)

perf-test-perf-responder: $(TARGET_PERF_MONT_ASM)
	$(PERF) ./$(TARGET_PERF_MONT_ASM) responder 20000

perf-test-perf-initiator: $(TARGET_PERF_MONT_ASM)
	$(PERF) ./$(TARGET_PERF_MONT_ASM) initiator 20000

perf-test-perf-scalar: $(TARGET_PERF_MONT_ASM)
	$(PERF) ./$(TARGET_PERF_MONT_ASM) scalar 50000

perf-test-perf-genra: $(TARGET_PERF_MONT_ASM)
	$(PERF) ./$(TARGET_PERF_MONT_ASM) genra 50000

# -------------------------------------------------
# Helpers
# -------------------------------------------------
print-flags:
	@echo "CC=$(CC)"
	@echo "CPU_TUNE=$(CPU_TUNE)"
	@echo "CPPFLAGS=$(CPPFLAGS)"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "ASFLAGS=$(ASFLAGS)"
	@echo "LDFLAGS=$(LDFLAGS)"

clean:
	rm -f \
	    $(TARGET_PRIME_C) \
	    $(TARGET_PRIME_FP_ASM) \
	    $(TARGET_MONT_C) \
	    $(TARGET_MONT_ASM) \
	    $(TARGET_PORNIN_C) \
	    $(TARGET_PORNIN_ASM) \
	    $(TARGET_PORNIN_INV_C) \
	    $(TARGET_PORNIN_INV_ASM) \
	    $(TARGET_PORNIN_INV_PROF_C) \
	    $(TARGET_PORNIN_INV_PROF_ASM) \
	    $(TARGET_INV_PERF) \
	    $(TARGET_INV_PERF_MONT_C) \
	    $(TARGET_INV_PERF_PORNIN_C) \
	    $(TARGET_PERF_MONT_ASM) \
	    $(TARGET_SECURE_MESSAGE_DEMO) \
	    *.o \
	    *.gcda *.gcno \
	    default.profraw
	rm -rf $(SECURE_BUILD_DIR)

.PHONY: all clean print-flags \
        run-prime-c run-prime-asm \
        run-mont-c run-mont-asm \
        run-pornin-c run-pornin-asm \
        run-pornin-inv-c run-pornin-inv-asm \
        run-pornin-inv-prof-c run-pornin-inv-prof-asm \
        run-test-perf run-secure-message-demo run-inv-perf \
        run-inv-perf-mont-c run-inv-perf-pornin-c \
        bench-prime bench-mont bench-all \
        perf-prime-c perf-mont-c perf-pornin-c \
        perf-prime-asm perf-mont-asm perf-pornin-asm perf-inv \
        perf-inv-mont-c perf-inv-pornin-c perf-all \
        perf-test-perf-responder perf-test-perf-initiator \
        perf-test-perf-scalar perf-test-perf-genra
