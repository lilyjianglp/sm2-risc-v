# =================================================================
# SM2 RISC-V 64 Optimization Project Makefile
# Muse Pi Pro / SpacemiT K1(X60) board-only tuned version
#
# Build matrix:
#   baseline-c    : Prime route + fp C + fn C, no asm
#   baseline-asm  : Prime route + fp ASM + fn ASM
#   mont-c        : Montgomery route + fp C + fn C, no asm
#   mont-asm      : Montgomery route + fp ASM + fn ASM
#
# Correctness matrix:
#   test_correct_prime_c       : baseline-c correctness
#   test_correct_prime_fp_asm  : baseline-asm correctness
#   test_correct_mont_c        : mont-c correctness
#   test_correct_mont_asm      : mont-asm correctness
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
# -----------------------------------------------------------------
CPU_TUNE   := spacemit-x60
CPU_CFLAGS := -mcpu=spacemit-x60

# -----------------------------------------------------------------
# Optimization / common flags
# -----------------------------------------------------------------
BASE_CFLAGS = -Wall -Wextra -Wshadow -Wno-unused-parameter -std=gnu11
OPT_CFLAGS  = -O3 -fomit-frame-pointer -fno-plt -falign-functions=32
SAFE_CFLAGS = -fno-strict-aliasing

# 汇编文件走预处理，方便 #ifdef / 宏展开
ASM_CPPFLAGS = -x assembler-with-cpp

ifeq ($(USE_LTO),1)
    LTO_FLAGS = -flto
else
    LTO_FLAGS =
endif

ifeq ($(USE_PGO),gen)
    PGO_FLAGS = -fprofile-generate
else ifeq ($(USE_PGO),use)
    PGO_FLAGS = -fprofile-use -fprofile-correction
else
    PGO_FLAGS =
endif

CPPFLAGS =

CFLAGS = \
    $(BASE_CFLAGS) \
    $(OPT_CFLAGS) \
    $(SAFE_CFLAGS) \
    $(CPU_CFLAGS) \
    $(LTO_FLAGS) \
    $(PGO_FLAGS)

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
# Test entry files
# -------------------------------------------------
TEST_SRC    = test.c
CORRECT_SRC = test_correct.c

# -------------------------------------------------
# fn layer variants
# -------------------------------------------------

# fn 纯 C 版本
# 用于 baseline-c / mont-c
FN_C_ONLY_C    = fn_Montgomery.c
FN_C_ONLY_DEFS = -DNO_TEST_FN_INTERNALS

# fn ASM 版本
# 用于 baseline-asm / mont-asm
FN_ASM_C       = fn_Montgomery.c
FN_ASM_S       = asm/fn_mont.S
FN_ASM_DEFS    = -DHAVE_FN_MONT_CORE -DUSE_ASM_FN_MONT -DNO_TEST_FN_INTERNALS

# -------------------------------------------------
# Prime / 普通模 p 路线
# -------------------------------------------------
FP_PRIME_C     = fp.c
CURVE_PRIME_C  = sm2_curve.c
SCALAR_PRIME_C = sm2_scalar.c
KEX_PRIME_C    = sm2_kex.c

PRIME_BASE_SRCS = \
    $(FP_PRIME_C) \
    $(CURVE_PRIME_C) \
    $(SCALAR_PRIME_C) \
    $(KEX_PRIME_C) \
    $(TEST_SRC)

PRIME_CORRECT_SRCS = \
    $(FP_PRIME_C) \
    $(CURVE_PRIME_C) \
    $(SCALAR_PRIME_C) \
    $(KEX_PRIME_C)

# -------------------------------------------------
# Montgomery 模 p 路线
# -------------------------------------------------
FP_MONT_C      = fp_montgomery.c
CURVE_MONT_C   = sm2_curve_mont.c
SCALAR_MONT_C  = sm2_scalar_mont.c
KEX_MONT_C     = sm2_kex_mont.c

MONT_BASE_SRCS = \
    $(FP_MONT_C) \
    $(CURVE_MONT_C) \
    $(SCALAR_MONT_C) \
    $(KEX_MONT_C) \
    $(TEST_SRC)

MONT_CORRECT_SRCS = \
    $(FP_MONT_C) \
    $(CURVE_MONT_C) \
    $(SCALAR_MONT_C) \
    $(KEX_MONT_C)

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

TARGET_PERF_MONT_ASM = test_perf

TARGET_CORRECT_PRIME_C      = test_correct_prime_c
TARGET_CORRECT_PRIME_FP_ASM = test_correct_prime_fp_asm
TARGET_CORRECT_MONT_C       = test_correct_mont_c
TARGET_CORRECT_MONT_ASM     = test_correct_mont_asm

all: $(TARGET_PRIME_C) $(TARGET_PRIME_FP_ASM) $(TARGET_MONT_C) $(TARGET_MONT_ASM)

# -------------------------------------------------
# baseline-c: Prime + fp C + fn C
# 不编译任何 .S 文件
# -------------------------------------------------
$(TARGET_PRIME_C): $(PRIME_BASE_SRCS) $(FN_C_ONLY_C)
	@echo ">>> Building baseline-c: Prime / fp C / fn C ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      $(FN_C_ONLY_DEFS) \
	      $(PRIME_BASE_SRCS) $(FN_C_ONLY_C) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# baseline-asm: Prime + fp ASM + fn ASM
# -------------------------------------------------
$(TARGET_PRIME_FP_ASM): $(PRIME_BASE_SRCS) $(FN_ASM_C) $(ASM_FP_PRIME_SRCS) $(FN_ASM_S)
	@echo ">>> Building baseline-asm: Prime / fp ASM / fn ASM ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      $(FN_ASM_DEFS) \
	      $(ASM_FP_PRIME_DEFS) \
	      $(PRIME_BASE_SRCS) $(FN_ASM_C) \
	      $(ASFLAGS) $(ASM_FP_PRIME_SRCS) $(FN_ASM_S) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# mont-c: Montgomery + fp C + fn C
# 不编译任何 .S 文件
# -------------------------------------------------
$(TARGET_MONT_C): $(MONT_BASE_SRCS) $(FN_C_ONLY_C)
	@echo ">>> Building mont-c: Montgomery / fp C / fn C ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT \
	      $(FN_C_ONLY_DEFS) \
	      $(MONT_BASE_SRCS) $(FN_C_ONLY_C) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# mont-asm: Montgomery + fp ASM + fn ASM
# -------------------------------------------------
$(TARGET_MONT_ASM): $(MONT_BASE_SRCS) $(FN_ASM_C) $(ASM_FP_MONT_SRCS) $(FN_ASM_S)
	@echo ">>> Building mont-asm: Montgomery / fp ASM / fn ASM ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT \
	      $(FN_ASM_DEFS) \
	      $(ASM_FP_MONT_DEFS) \
	      $(MONT_BASE_SRCS) $(FN_ASM_C) \
	      $(ASFLAGS) $(ASM_FP_MONT_SRCS) $(FN_ASM_S) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# test_perf: 默认使用 mont-asm 路线
# -------------------------------------------------
MONT_PERF_SRCS = $(filter-out $(TEST_SRC),$(MONT_BASE_SRCS))

$(TARGET_PERF_MONT_ASM): test_perf.c $(MONT_PERF_SRCS) $(FN_ASM_C) $(ASM_FP_MONT_SRCS) $(FN_ASM_S)
	@echo ">>> Building test_perf: Montgomery / fp ASM / fn ASM ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT \
	      $(FN_ASM_DEFS) \
	      $(ASM_FP_MONT_DEFS) \
	      test_perf.c \
	      $(MONT_PERF_SRCS) $(FN_ASM_C) \
	      $(ASFLAGS) $(ASM_FP_MONT_SRCS) $(FN_ASM_S) \
	      $(LDFLAGS) -o $@

# -------------------------------------------------
# Correctness tests
# 注意：
# 这里不要把 Prime 和 Montgomery 源码混在一个程序里链接，
# 否则会出现 multiple definition。
# 每个正确性程序只链接自己对应路线的源码。
# -------------------------------------------------

$(TARGET_CORRECT_PRIME_C): $(CORRECT_SRC) $(PRIME_CORRECT_SRCS) $(FN_C_ONLY_C)
	@echo ">>> Building correctness: baseline-c / Prime / fp C / fn C ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      $(FN_C_ONLY_DEFS) \
	      $(CORRECT_SRC) \
	      $(PRIME_CORRECT_SRCS) $(FN_C_ONLY_C) \
	      $(LDFLAGS) -o $@

$(TARGET_CORRECT_PRIME_FP_ASM): $(CORRECT_SRC) $(PRIME_CORRECT_SRCS) $(FN_ASM_C) $(ASM_FP_PRIME_SRCS) $(FN_ASM_S)
	@echo ">>> Building correctness: baseline-asm / Prime / fp ASM / fn ASM ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      $(FN_ASM_DEFS) \
	      $(ASM_FP_PRIME_DEFS) \
	      $(CORRECT_SRC) \
	      $(PRIME_CORRECT_SRCS) $(FN_ASM_C) \
	      $(ASFLAGS) $(ASM_FP_PRIME_SRCS) $(FN_ASM_S) \
	      $(LDFLAGS) -o $@

$(TARGET_CORRECT_MONT_C): $(CORRECT_SRC) $(MONT_CORRECT_SRCS) $(FN_C_ONLY_C)
	@echo ">>> Building correctness: mont-c / Montgomery / fp C / fn C ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT \
	      $(FN_C_ONLY_DEFS) \
	      $(CORRECT_SRC) \
	      $(MONT_CORRECT_SRCS) $(FN_C_ONLY_C) \
	      $(LDFLAGS) -o $@

$(TARGET_CORRECT_MONT_ASM): $(CORRECT_SRC) $(MONT_CORRECT_SRCS) $(FN_ASM_C) $(ASM_FP_MONT_SRCS) $(FN_ASM_S)
	@echo ">>> Building correctness: mont-asm / Montgomery / fp ASM / fn ASM ..."
	$(CC) $(CPPFLAGS) $(CFLAGS) \
	      -DUSE_FP_MONT \
	      $(FN_ASM_DEFS) \
	      $(ASM_FP_MONT_DEFS) \
	      $(CORRECT_SRC) \
	      $(MONT_CORRECT_SRCS) $(FN_ASM_C) \
	      $(ASFLAGS) $(ASM_FP_MONT_SRCS) $(FN_ASM_S) \
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

run-correct-prime-c: $(TARGET_CORRECT_PRIME_C)
	./$(TARGET_CORRECT_PRIME_C)

run-correct-prime-asm: $(TARGET_CORRECT_PRIME_FP_ASM)
	./$(TARGET_CORRECT_PRIME_FP_ASM)

run-correct-mont-c: $(TARGET_CORRECT_MONT_C)
	./$(TARGET_CORRECT_MONT_C)

run-correct-mont-asm: $(TARGET_CORRECT_MONT_ASM)
	./$(TARGET_CORRECT_MONT_ASM)

# -------------------------------------------------
# User-friendly aliases
# -------------------------------------------------
baseline-c: run-prime-c

baseline-asm: run-prime-asm

mont-c: run-mont-c

mont-asm: run-mont-asm

demo: mont-asm

test: correctness

bench: bench-all

correctness: $(TARGET_CORRECT_PRIME_C) $(TARGET_CORRECT_PRIME_FP_ASM) $(TARGET_CORRECT_MONT_C) $(TARGET_CORRECT_MONT_ASM)
	@echo "========================================================"
	@echo "=== correctness: baseline-c / Prime / fp C / fn C ==="
	./$(TARGET_CORRECT_PRIME_C)
	@echo ""
	@echo "=== correctness: baseline-asm / Prime / fp ASM / fn ASM ==="
	./$(TARGET_CORRECT_PRIME_FP_ASM)
	@echo ""
	@echo "=== correctness: mont-c / Montgomery / fp C / fn C ==="
	./$(TARGET_CORRECT_MONT_C)
	@echo ""
	@echo "=== correctness: mont-asm / Montgomery / fp ASM / fn ASM ==="
	./$(TARGET_CORRECT_MONT_ASM)
	@echo "========================================================"
	@echo "[PASS] correctness workflow completed"

# -------------------------------------------------
# Bench
# -------------------------------------------------
bench-prime: $(TARGET_PRIME_C) $(TARGET_PRIME_FP_ASM)
	@echo "========================================================"
	@echo "CPU_TUNE=$(CPU_TUNE)"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "ASFLAGS=$(ASFLAGS)"
	@echo ""
	@echo "=== baseline-c: Prime / fp C / fn C ==="
	./$(TARGET_PRIME_C)
	@echo ""
	@echo "=== baseline-asm: Prime / fp ASM / fn ASM ==="
	./$(TARGET_PRIME_FP_ASM)
	@echo "========================================================"

bench-mont: $(TARGET_MONT_C) $(TARGET_MONT_ASM)
	@echo "========================================================"
	@echo "CPU_TUNE=$(CPU_TUNE)"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "ASFLAGS=$(ASFLAGS)"
	@echo ""
	@echo "=== mont-c: Montgomery / fp C / fn C ==="
	./$(TARGET_MONT_C)
	@echo ""
	@echo "=== mont-asm: Montgomery / fp ASM / fn ASM ==="
	./$(TARGET_MONT_ASM)
	@echo "========================================================"

bench-all: $(TARGET_PRIME_C) $(TARGET_PRIME_FP_ASM) $(TARGET_MONT_C) $(TARGET_MONT_ASM)
	@echo "========================================================"
	@echo "CPU_TUNE=$(CPU_TUNE)"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "ASFLAGS=$(ASFLAGS)"
	@echo ""
	@echo "=== baseline-c: Prime / fp C / fn C ==="
	./$(TARGET_PRIME_C)
	@echo ""
	@echo "=== baseline-asm: Prime / fp ASM / fn ASM ==="
	./$(TARGET_PRIME_FP_ASM)
	@echo ""
	@echo "=== mont-c: Montgomery / fp C / fn C ==="
	./$(TARGET_MONT_C)
	@echo ""
	@echo "=== mont-asm: Montgomery / fp ASM / fn ASM ==="
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
	@echo "=== baseline-c: Prime / fp C / fn C ==="
	$(PERF) ./$(TARGET_PRIME_C)
	@echo ""
	@echo "=== baseline-asm: Prime / fp ASM / fn ASM ==="
	$(PERF) ./$(TARGET_PRIME_FP_ASM)
	@echo ""
	@echo "=== mont-c: Montgomery / fp C / fn C ==="
	$(PERF) ./$(TARGET_MONT_C)
	@echo ""
	@echo "=== mont-asm: Montgomery / fp ASM / fn ASM ==="
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
	@echo ""
	@echo "Build matrix:"
	@echo "  baseline-c   = Prime + fp C + fn C"
	@echo "  baseline-asm = Prime + fp ASM + fn ASM"
	@echo "  mont-c       = Montgomery + fp C + fn C"
	@echo "  mont-asm     = Montgomery + fp ASM + fn ASM"
	@echo ""
	@echo "Correctness:"
	@echo "  make correctness = run SM2 correctness tests for all four variants"

clean:
	rm -f \
	    $(TARGET_PRIME_C) \
	    $(TARGET_PRIME_FP_ASM) \
	    $(TARGET_MONT_C) \
	    $(TARGET_MONT_ASM) \
	    $(TARGET_PERF_MONT_ASM) \
	    $(TARGET_CORRECT_PRIME_C) \
	    $(TARGET_CORRECT_PRIME_FP_ASM) \
	    $(TARGET_CORRECT_MONT_C) \
	    $(TARGET_CORRECT_MONT_ASM) \
	    *.o \
	    *.gcda *.gcno \
	    default.profraw

.PHONY: all clean print-flags \
        run-prime-c run-prime-asm run-mont-c run-mont-asm run-test-perf \
        run-correct-prime-c run-correct-prime-asm run-correct-mont-c run-correct-mont-asm \
        baseline-c baseline-asm mont-c mont-asm demo test bench correctness \
        bench-prime bench-mont bench-all \
        perf-prime-c perf-mont-c perf-prime-asm perf-mont-asm perf-all \
        perf-test-perf-responder perf-test-perf-initiator \
        perf-test-perf-scalar perf-test-perf-genra