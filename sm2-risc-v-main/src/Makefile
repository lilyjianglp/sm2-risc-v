# =================================================================
# SM2 RISC-V 64 Optimization Project Makefile
# Muse Pi Pro / SpacemiT K1(X60) board-only tuned version
#
# fn 层统一固定为：fn_Montgomery.c + asm/fn_mont.S
# 这样可以保证 fn_mont_mul_ct 等符号由 asm/fn_mont.S 提供
# =================================================================

CC      ?= gcc
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
# fn 层：四个 target 全部固定
# 恢复 asm/fn_mont.S，以提供 fn_mont_mul_ct 等符号
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
# Montgomery 模 p 路线
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
TARGET_PRIME_C      = test_prime_c
TARGET_PRIME_FP_ASM = test_prime_fp_asm
TARGET_MONT_C       = test_mont_c
TARGET_MONT_ASM     = test_mont_asm

all: $(TARGET_PRIME_C) $(TARGET_PRIME_FP_ASM) $(TARGET_MONT_C) $(TARGET_MONT_ASM)

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
# Montgomery 纯 C（fn 固定为 C+asm）
# -------------------------------------------------
$(TARGET_MONT_C): $(MONT_SRCS) $(FN_FIXED_ASM)
	@echo ">>> Building Montgomery / Pure C (fn fixed to C+asm) ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT $(FN_FIXED_DEFS) \
	      $(MONT_SRCS) $(FN_FIXED_ASM) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# Montgomery ASM（fn 固定为 C+asm）
# -------------------------------------------------
$(TARGET_MONT_ASM): $(MONT_SRCS) $(ASM_FP_MONT_SRCS) $(FN_FIXED_ASM)
	@echo ">>> Building Montgomery / fp ASM (fn fixed to C+asm) ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT $(FN_FIXED_DEFS) \
	      $(ASM_FP_MONT_DEFS) \
	      $(MONT_SRCS) \
	      $(ASFLAGS) $(ASM_FP_MONT_SRCS) $(FN_FIXED_ASM) \
	      $(LDFLAGS) -o $@


TARGET_PERF_MONT_ASM = test_perf
MONT_PERF_SRCS = $(filter-out $(TEST_SRC),$(MONT_SRCS))

$(TARGET_PERF_MONT_ASM): test_perf.c $(MONT_PERF_SRCS) $(ASM_FP_MONT_SRCS) $(FN_FIXED_ASM)
	@echo ">>> Building test_perf (Montgomery / fp ASM, fn fixed to C+asm) ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT $(FN_FIXED_DEFS) \
	      $(ASM_FP_MONT_DEFS) \
	      test_perf.c \
	      $(MONT_PERF_SRCS) \
	      $(ASFLAGS) $(ASM_FP_MONT_SRCS) $(FN_FIXED_ASM) \
	      $(LDFLAGS) -o $@


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
run-test-perf: $(TARGET_PERF_MONT_ASM)
	./$(TARGET_PERF_MONT_ASM)
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

bench-mont: $(TARGET_MONT_C) $(TARGET_MONT_ASM)
	@echo "========================================================"
	@echo "CPU_TUNE=$(CPU_TUNE)  fn=C+asm fixed"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "ASFLAGS=$(ASFLAGS)"
	@echo "=== Montgomery / Pure C ==="
	./$(TARGET_MONT_C)
	@echo ""
	@echo "=== Montgomery / fp ASM ==="
	./$(TARGET_MONT_ASM)
	@echo "========================================================"

bench-all: $(TARGET_PRIME_C) $(TARGET_PRIME_FP_ASM) $(TARGET_MONT_C) $(TARGET_MONT_ASM)
	@echo "========================================================"
	@echo "CPU_TUNE=$(CPU_TUNE)  fn=C+asm fixed"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "ASFLAGS=$(ASFLAGS)"
	@echo "=== Prime / Pure C ==="
	./$(TARGET_PRIME_C)
	@echo ""
	@echo "=== Montgomery / Pure C ==="
	./$(TARGET_MONT_C)
	@echo ""
	@echo "=== Prime / fp ASM ==="
	./$(TARGET_PRIME_FP_ASM)
	@echo ""
	@echo "=== Montgomery / fp ASM ==="
	./$(TARGET_MONT_ASM)
	@echo "========================================================"

# -------------------------------------------------
# Perf
# -------------------------------------------------
PERF = perf stat -e cycles,instructions

perf-prime-c: $(TARGET_PRIME_C)
	$(PERF) ./$(TARGET_PRIME_C)

perf-mont-c: $(TARGET_MONT_C)
	$(PERF) ./$(TARGET_MONT_C)

perf-prime-asm: $(TARGET_PRIME_FP_ASM)
	$(PERF) ./$(TARGET_PRIME_FP_ASM)

perf-mont-asm: $(TARGET_MONT_ASM)
	$(PERF) ./$(TARGET_MONT_ASM)

perf-all: $(TARGET_PRIME_C) $(TARGET_PRIME_FP_ASM) $(TARGET_MONT_C) $(TARGET_MONT_ASM)
	@echo "=== Prime / Pure C ==="
	$(PERF) ./$(TARGET_PRIME_C)
	@echo ""
	@echo "=== Montgomery / Pure C ==="
	$(PERF) ./$(TARGET_MONT_C)
	@echo ""
	@echo "=== Prime / fp ASM ==="
	$(PERF) ./$(TARGET_PRIME_FP_ASM)
	@echo ""
	@echo "=== Montgomery / fp ASM ==="
	$(PERF) ./$(TARGET_MONT_ASM)
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
		$(TARGET_PERF_MONT_ASM) \
	    *.o \
	    *.gcda *.gcno \
	    default.profraw

.PHONY: all clean print-flags \
        run-prime-c run-prime-asm run-mont-c run-mont-asm run-test-perf \
        bench-prime bench-mont bench-all \
        perf-prime-c perf-mont-c perf-prime-asm perf-mont-asm perf-all \
        perf-test-perf-responder perf-test-perf-initiator \
        perf-test-perf-scalar perf-test-perf-genra