CC      := gcc
OBJCOPY := objcopy

SM2_SRC   := /home/sm2/sm2-risc-v-main/sm2-risc-v-main/src
GMSSL_INC := /home/sm2/GmSSL/include
GMSSL_LIB := /home/sm2/GmSSL/build/bin

# 关键：关闭 stack protector，避免 objcopy 前缀化 __stack_chk_* 导致链接失败
CFLAGS_COMMON := -O3 -fno-stack-protector -I$(SM2_SRC) -I$(GMSSL_INC)
LDFLAGS_COMMON := -L$(GMSSL_LIB) -lgmssl

BUILD_DIR    := build
C_DIR        := $(BUILD_DIR)/c
ASM_DIR      := $(BUILD_DIR)/asm
C_PREF_DIR   := $(BUILD_DIR)/c_pref
ASM_PREF_DIR := $(BUILD_DIR)/asm_pref

# =========================================================
# 未来 bench_all.c 扩展时常用到的 core C 源文件
# 这几份一次性都编好，后面 bench_all.c 基本不用再改 Makefile
# =========================================================
CORE_C_SRCS := \
	$(SM2_SRC)/fp_montgomery.c \
	$(SM2_SRC)/fn_Montgomery.c \
	$(SM2_SRC)/sm2_curve_mont.c \
	$(SM2_SRC)/sm2_scalar_mont.c

# =========================================================
# Fp 相关 ASM 源文件
# 以后加 fp add/sub/sqr/to/from 基本都不用改这里
# =========================================================
FP_ASM_SRCS := \
	$(SM2_SRC)/asm/fp_mont_add.S \
	$(SM2_SRC)/asm/fp_mont_sub.S \
	$(SM2_SRC)/asm/fp_mont_neg.S \
	$(SM2_SRC)/asm/fp_mont_mul.S \
	$(SM2_SRC)/asm/fp_mont_sqr.S \
	$(SM2_SRC)/asm/fp_mont_to_from.S

FN_ASM_SRCS := \
	$(SM2_SRC)/asm/fn_mont.S

# =========================================================
# 原始对象文件
# =========================================================
C_CORE_OBJS := \
	$(C_DIR)/fp_montgomery.o \
	$(C_DIR)/fn_Montgomery.o \
	$(C_DIR)/sm2_curve_mont.o \
	$(C_DIR)/sm2_scalar_mont.o

ASM_CORE_OBJS := \
	$(ASM_DIR)/fp_montgomery.o \
	$(ASM_DIR)/fn_Montgomery.o \
	$(ASM_DIR)/sm2_curve_mont.o \
	$(ASM_DIR)/sm2_scalar_mont.o

ASM_FP_OBJS := \
	$(ASM_DIR)/fp_mont_add.o \
	$(ASM_DIR)/fp_mont_sub.o \
	$(ASM_DIR)/fp_mont_neg.o \
	$(ASM_DIR)/fp_mont_mul.o \
	$(ASM_DIR)/fp_mont_sqr.o \
	$(ASM_DIR)/fp_mont_to_from.o

ASM_FN_OBJS := \
	$(ASM_DIR)/fn_mont.o

# =========================================================
# 前缀化后的对象文件
# =========================================================
C_CORE_PREF_OBJS := \
	$(C_PREF_DIR)/fp_montgomery.o \
	$(C_PREF_DIR)/fn_Montgomery.o \
	$(C_PREF_DIR)/sm2_curve_mont.o \
	$(C_PREF_DIR)/sm2_scalar_mont.o

ASM_CORE_PREF_OBJS := \
	$(ASM_PREF_DIR)/fp_montgomery.o \
	$(ASM_PREF_DIR)/fn_Montgomery.o \
	$(ASM_PREF_DIR)/sm2_curve_mont.o \
	$(ASM_PREF_DIR)/sm2_scalar_mont.o

ASM_FP_PREF_OBJS := \
	$(ASM_PREF_DIR)/fp_mont_add.o \
	$(ASM_PREF_DIR)/fp_mont_sub.o \
	$(ASM_PREF_DIR)/fp_mont_neg.o \
	$(ASM_PREF_DIR)/fp_mont_mul.o \
	$(ASM_PREF_DIR)/fp_mont_sqr.o \
	$(ASM_PREF_DIR)/fp_mont_to_from.o

ASM_FN_PREF_OBJS := \
	$(ASM_PREF_DIR)/fn_mont.o

ALL_PREF_OBJS := \
	$(C_CORE_PREF_OBJS) \
	$(ASM_CORE_PREF_OBJS) \
	$(ASM_FP_PREF_OBJS) \
	$(ASM_FN_PREF_OBJS)
.PHONY: all run perf clean dirs show

all: bench_all

run: bench_all
	env LD_LIBRARY_PATH=$(GMSSL_LIB) ./bench_all

perf: bench_all
	env LD_LIBRARY_PATH=$(GMSSL_LIB) perf stat -e cycles,instructions ./bench_all

show:
	@echo "SM2_SRC=$(SM2_SRC)"
	@echo "GMSSL_INC=$(GMSSL_INC)"
	@echo "GMSSL_LIB=$(GMSSL_LIB)"
	@echo "CFLAGS_COMMON=$(CFLAGS_COMMON)"

dirs:
	mkdir -p $(C_DIR) $(ASM_DIR) $(C_PREF_DIR) $(ASM_PREF_DIR)

# =========================================================
# 1) 编原始 C 版本对象
# =========================================================
$(C_DIR)/fp_montgomery.o: $(SM2_SRC)/fp_montgomery.c | dirs
	$(CC) $(CFLAGS_COMMON) -DBENCH_USE_ASM=0 -c $< -o $@

$(C_DIR)/fn_Montgomery.o: $(SM2_SRC)/fn_Montgomery.c | dirs
	$(CC) $(CFLAGS_COMMON) -DBENCH_USE_ASM=0 -c $< -o $@

$(C_DIR)/sm2_curve_mont.o: $(SM2_SRC)/sm2_curve_mont.c | dirs
	$(CC) $(CFLAGS_COMMON) -DBENCH_USE_ASM=0 -c $< -o $@

$(C_DIR)/sm2_scalar_mont.o: $(SM2_SRC)/sm2_scalar_mont.c | dirs
	$(CC) $(CFLAGS_COMMON) -DBENCH_USE_ASM=0 -c $< -o $@

# =========================================================
# 2) 编原始 ASM 版本对象
# =========================================================
$(ASM_DIR)/fp_montgomery.o: $(SM2_SRC)/fp_montgomery.c | dirs
	$(CC) $(CFLAGS_COMMON) \
		-DBENCH_USE_ASM=1 \
		-DUSE_RV64_ASM_ADD \
		-DUSE_RV64_ASM_SUB \
		-DUSE_RV64_ASM_NEG \
		-DUSE_RV64_ASM_MUL \
		-DUSE_RV64_ASM_SQR \
		-DUSE_RV64_ASM_TO_FROM \
		-c $< -o $@

$(ASM_DIR)/fn_Montgomery.o: $(SM2_SRC)/fn_Montgomery.c | dirs
	$(CC) $(CFLAGS_COMMON) \
		-DBENCH_USE_ASM=1 \
		-DUSE_RV64_ASM_ADD \
		-DUSE_RV64_ASM_SUB \
		-DUSE_RV64_ASM_NEG \
		-DUSE_RV64_ASM_MUL \
		-DUSE_RV64_ASM_SQR \
		-DUSE_RV64_ASM_TO_FROM \
		-c $< -o $@

$(ASM_DIR)/sm2_curve_mont.o: $(SM2_SRC)/sm2_curve_mont.c | dirs
	$(CC) $(CFLAGS_COMMON) \
		-DBENCH_USE_ASM=1 \
		-DUSE_RV64_ASM_ADD \
		-DUSE_RV64_ASM_SUB \
		-DUSE_RV64_ASM_NEG \
		-DUSE_RV64_ASM_MUL \
		-DUSE_RV64_ASM_SQR \
		-DUSE_RV64_ASM_TO_FROM \
		-c $< -o $@

$(ASM_DIR)/sm2_scalar_mont.o: $(SM2_SRC)/sm2_scalar_mont.c | dirs
	$(CC) $(CFLAGS_COMMON) \
		-DBENCH_USE_ASM=1 \
		-DUSE_RV64_ASM_ADD \
		-DUSE_RV64_ASM_SUB \
		-DUSE_RV64_ASM_NEG \
		-DUSE_RV64_ASM_MUL \
		-DUSE_RV64_ASM_SQR \
		-DUSE_RV64_ASM_TO_FROM \
		-c $< -o $@

# =========================================================
# 3) 编 ASM .S 对象
# =========================================================
$(ASM_DIR)/fp_mont_add.o: $(SM2_SRC)/asm/fp_mont_add.S | dirs
	$(CC) $(CFLAGS_COMMON) -c $< -o $@

$(ASM_DIR)/fp_mont_sub.o: $(SM2_SRC)/asm/fp_mont_sub.S | dirs
	$(CC) $(CFLAGS_COMMON) -c $< -o $@

$(ASM_DIR)/fp_mont_neg.o: $(SM2_SRC)/asm/fp_mont_neg.S | dirs
	$(CC) $(CFLAGS_COMMON) -c $< -o $@

$(ASM_DIR)/fp_mont_mul.o: $(SM2_SRC)/asm/fp_mont_mul.S | dirs
	$(CC) $(CFLAGS_COMMON) -c $< -o $@

$(ASM_DIR)/fp_mont_sqr.o: $(SM2_SRC)/asm/fp_mont_sqr.S | dirs
	$(CC) $(CFLAGS_COMMON) -c $< -o $@

$(ASM_DIR)/fp_mont_to_from.o: $(SM2_SRC)/asm/fp_mont_to_from.S | dirs
	$(CC) $(CFLAGS_COMMON) -c $< -o $@

$(ASM_DIR)/fn_mont.o: $(SM2_SRC)/asm/fn_mont.S | dirs
	$(CC) $(CFLAGS_COMMON) -c $< -o $@

# =========================================================
# 4) 给 C 版本对象统一加 c_ 前缀
# =========================================================
$(C_PREF_DIR)/fp_montgomery.o: $(C_DIR)/fp_montgomery.o | dirs
	$(OBJCOPY) --prefix-symbols=c_ $< $@

$(C_PREF_DIR)/fn_Montgomery.o: $(C_DIR)/fn_Montgomery.o | dirs
	$(OBJCOPY) --prefix-symbols=c_ $< $@

$(C_PREF_DIR)/sm2_curve_mont.o: $(C_DIR)/sm2_curve_mont.o | dirs
	$(OBJCOPY) --prefix-symbols=c_ $< $@

$(C_PREF_DIR)/sm2_scalar_mont.o: $(C_DIR)/sm2_scalar_mont.o | dirs
	$(OBJCOPY) --prefix-symbols=c_ $< $@


# =========================================================
# 5) 给 ASM 版本对象统一加 asm_ 前缀
# =========================================================
$(ASM_PREF_DIR)/fp_montgomery.o: $(ASM_DIR)/fp_montgomery.o | dirs
	$(OBJCOPY) --prefix-symbols=asm_ $< $@

$(ASM_PREF_DIR)/fn_Montgomery.o: $(ASM_DIR)/fn_Montgomery.o | dirs
	$(OBJCOPY) --prefix-symbols=asm_ $< $@

$(ASM_PREF_DIR)/sm2_curve_mont.o: $(ASM_DIR)/sm2_curve_mont.o | dirs
	$(OBJCOPY) --prefix-symbols=asm_ $< $@

$(ASM_PREF_DIR)/sm2_scalar_mont.o: $(ASM_DIR)/sm2_scalar_mont.o | dirs
	$(OBJCOPY) --prefix-symbols=asm_ $< $@

$(ASM_PREF_DIR)/fp_mont_add.o: $(ASM_DIR)/fp_mont_add.o | dirs
	$(OBJCOPY) --prefix-symbols=asm_ $< $@

$(ASM_PREF_DIR)/fp_mont_sub.o: $(ASM_DIR)/fp_mont_sub.o | dirs
	$(OBJCOPY) --prefix-symbols=asm_ $< $@

$(ASM_PREF_DIR)/fp_mont_neg.o: $(ASM_DIR)/fp_mont_neg.o | dirs
	$(OBJCOPY) --prefix-symbols=asm_ $< $@

$(ASM_PREF_DIR)/fp_mont_mul.o: $(ASM_DIR)/fp_mont_mul.o | dirs
	$(OBJCOPY) --prefix-symbols=asm_ $< $@

$(ASM_PREF_DIR)/fp_mont_sqr.o: $(ASM_DIR)/fp_mont_sqr.o | dirs
	$(OBJCOPY) --prefix-symbols=asm_ $< $@

$(ASM_PREF_DIR)/fp_mont_to_from.o: $(ASM_DIR)/fp_mont_to_from.o | dirs
	$(OBJCOPY) --prefix-symbols=asm_ $< $@

$(ASM_PREF_DIR)/fn_mont.o: $(ASM_DIR)/fn_mont.o | dirs
	$(OBJCOPY) --prefix-symbols=asm_ $< $@

# =========================================================
# 6) 最终链接统一 bench 程序
# =========================================================
bench_all: bench_all.c $(ALL_PREF_OBJS)
	$(CC) $(CFLAGS_COMMON) \
		bench_all.c \
		$(ALL_PREF_OBJS) \
		$(LDFLAGS_COMMON) \
		-o $@

clean:
	rm -rf $(BUILD_DIR) bench_all