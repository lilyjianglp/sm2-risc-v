# =================================================================
# SM2 RISC-V 64 Optimization Project Makefile
# =================================================================

# 编译器设置
CC      = riscv64-linux-gnu-gcc
CFLAGS  = -O2 -static -Wall

# 模拟器设置 (用于运行测试)
SIM     = qemu-riscv64

# -----------------------------------------------------------------
# 文件列表
# -----------------------------------------------------------------

# 通用 C 源码 (基础库 + 上层逻辑 + 测试入口)
C_SRCS  = fp.c \
          fn.c \
          sm2_curve.c \
          sm2_scalar.c \
          sm2_kex.c \
          test.c

# 汇编源码 (在 asm/ 目录下)
ASM_SRCS = asm/fp_add.S \
           asm/fp_sub.S \
           asm/fp_neg.S \
           asm/fp_mul.S \
           asm/fp_sqr.S \
           asm/fp_reduce.S

# 汇编优化宏定义 (开启所有汇编加速)
ASM_FLAGS = -DUSE_ASM_FP_ADD \
            -DUSE_ASM_FP_SUB \
            -DUSE_ASM_FP_NEG \
            -DUSE_ASM_FP_MUL \
            -DUSE_ASM_FP_SQR \
            -DUSE_ASM_REDUCE

# -----------------------------------------------------------------
# 编译目标
# -----------------------------------------------------------------

# 输出文件名
TARGET_C  = test_sm2_kex_c
TARGET_RV = test_sm2_kex_rv

# 默认目标：编译两个版本
all: $(TARGET_C) $(TARGET_RV)

# 1. 编译纯 C 版本 (Baseline)
# 不包含 ASM_SRCS，也不包含 ASM_FLAGS
$(TARGET_C): $(C_SRCS)
	@echo "正在编译纯 C 版本: $@"
	$(CC) $(CFLAGS) $(C_SRCS) -o $@

# 2. 编译汇编优化版本 (Optimized)
# 包含 ASM_SRCS，并定义 ASM_FLAGS 宏
$(TARGET_RV): $(C_SRCS) $(ASM_SRCS)
	@echo "正在编译汇编优化版本: $@"
	$(CC) $(CFLAGS) $(ASM_FLAGS) $(C_SRCS) $(ASM_SRCS) -o $@

# -----------------------------------------------------------------
# 辅助命令
# -----------------------------------------------------------------

# 运行测试 (需要安装 qemu-riscv64)
run: all
	@echo "========================================================"
	@echo ">>> 运行纯 C 版本测试 (Baseline)..."
	$(SIM) ./$(TARGET_C)
	@echo ""
	@echo "========================================================"
	@echo ">>> 运行汇编优化版本测试 (Optimized)..."
	$(SIM) ./$(TARGET_RV)
	@echo "========================================================"

# 清理生成的文件
clean:
	rm -f $(TARGET_C) $(TARGET_RV) *.o

.PHONY: all run clean