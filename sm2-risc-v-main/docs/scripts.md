# 一键脚本说明

本项目新增 `scripts/` 目录，用于提供一键复现脚本，方便快速运行正确性测试、性能测试和代码行数统计。

当前脚本包括：

```text
scripts/
├── run_correctness.sh
├── run_perf.sh
└── count_loc.sh
```

---

## 1. run_correctness.sh

### 作用

`run_correctness.sh` 用于一键运行 SM2 正确性测试流程。

该脚本会进入 `src/` 目录，清理旧编译产物，然后执行：

```bash
make clean
make correctness
```

### 运行命令

在项目根目录执行：

```bash
bash scripts/run_correctness.sh
```

### 对应测试内容

该脚本会依次构建并运行四个正确性测试程序：

```text
test_correct_prime_c
test_correct_prime_fp_asm
test_correct_mont_c
test_correct_mont_asm
```

对应四个构建版本：

| 测试程序 | 对应版本 |
|---|---|
| `test_correct_prime_c` | `baseline-c` |
| `test_correct_prime_fp_asm` | `baseline-asm` |
| `test_correct_mont_c` | `mont-c` |
| `test_correct_mont_asm` | `mont-asm` |

每个版本都会执行：

- 固定 SM2 KEX 测试；
- Alice / Bob 共享密钥一致性检查；
- S1 / S2 确认值检查；
- 10 轮随机 SM2 KEX 测试。

### 当前结果

四个版本均通过正确性测试，并且固定测试向量下输出完全一致：

```text
SharedKey = db0e9fdca7c9b8cde2edf89e6ec184ffeb3c34dbe50bc8ad72c11f7735f4108f
S1        = 53752ebb074279b5c909d88c4ae9bb611e5dc90dba3eb079ce1133d62ff489d3
S2        = 431f7b64ceaa65f8cadffcc5b6056530c2139b4f2fd26b2fb6d05e251c709b85
```

这说明 Montgomery 优化和 RISC-V ASM 优化没有破坏 SM2 密钥交换协议正确性。

---

## 2. run_perf.sh

### 作用

`run_perf.sh` 用于一键运行 Linux `perf` 性能测试流程。

该脚本会进入 `src/` 目录，开启板子上的 `perf` 用户访问权限，然后执行：

```bash
make clean
make all
make perf-all
```

### 运行命令

在项目根目录执行：

```bash
bash scripts/run_perf.sh
```

### 内部执行流程

脚本主要执行以下命令：

```bash
cd src
sudo sysctl -w kernel.perf_user_access=2
make clean
make all
make perf-all
```

### 输出内容

`make perf-all` 会依次对四个版本执行：

```bash
perf stat -e cycles,instructions ./test_prime_c
perf stat -e cycles,instructions ./test_prime_fp_asm
perf stat -e cycles,instructions ./test_mont_c
perf stat -e cycles,instructions ./test_mont_asm
```

对应四个构建版本：

| 可执行文件 | 对应版本 |
|---|---|
| `test_prime_c` | `baseline-c` |
| `test_prime_fp_asm` | `baseline-asm` |
| `test_mont_c` | `mont-c` |
| `test_mont_asm` | `mont-asm` |

输出指标包括：

- `cycles`；
- `instructions`；
- `IPC`；
- `elapsed time`。

### 当前性能结果

| 版本 | Cycles | Instructions | IPC | Time | 相对 baseline-c |
|---|---:|---:|---:|---:|---:|
| `baseline-c` | 19.40B | 31.92B | 1.64 | 12.13s | 1.00x |
| `baseline-asm` | 24.61B | 35.61B | 1.45 | 15.38s | 0.79x |
| `mont-c` | 11.12B | 14.32B | 1.29 | 6.95s | 1.75x |
| `mont-asm` | 9.55B | 12.80B | 1.34 | 5.97s | 2.03x |

结论：

- `mont-c` 相比 `baseline-c` 已经获得明显提升；
- `mont-asm` 是当前最快版本；
- `mont-asm` 相比 `baseline-c` 总运行时间从 12.13s 降低到 5.97s，整体约 2.03x 加速。

---

## 3. count_loc.sh

### 作用

`count_loc.sh` 用于统计项目源码、汇编代码、脚本和 Makefile 的代码行数。

该脚本主要用于展示项目代码规模，后续可以放入 README、项目报告或答辩材料中。

### 运行命令

在项目根目录执行：

```bash
bash scripts/count_loc.sh
```

### 统计范围

当前统计范围包括：

```text
src/
scripts/
```

统计文件类型包括：

```text
.c
.h
.S
.s
.sh
Makefile
```

### 输出内容

脚本会统计项目中主要源代码文件的行数，并输出每类文件或每个文件的代码规模。

典型用途包括：

- 展示项目总代码量；
- 区分 C 源码、头文件、汇编文件和脚本文件；
- 为 README、实验报告或答辩材料提供代码规模依据。

---

## 4. 脚本使用建议

建议在项目根目录下执行所有脚本：

```bash
bash scripts/run_correctness.sh
bash scripts/run_perf.sh
bash scripts/count_loc.sh
```

不要在 `scripts/` 目录内部直接执行脚本，否则可能因为相对路径不同导致找不到 `src/` 目录。

如果脚本没有执行权限，可以先添加权限：

```bash
chmod +x scripts/run_correctness.sh
chmod +x scripts/run_perf.sh
chmod +x scripts/count_loc.sh
```

之后也可以直接执行：

```bash
./scripts/run_correctness.sh
./scripts/run_perf.sh
./scripts/count_loc.sh
```

---

## 5. 总结

`scripts/` 目录提供了三个一键脚本：

| 脚本 | 作用 |
|---|---|
| `run_correctness.sh` | 一键运行 SM2 正确性测试 |
| `run_perf.sh` | 一键运行 Linux `perf` 性能测试 |
| `count_loc.sh` | 一键统计项目代码行数 |

通过这些脚本，可以快速复现项目的核心实验结果：

- 验证四个版本的 SM2 密钥交换正确性；
- 对比 `baseline-c`、`baseline-asm`、`mont-c`、`mont-asm` 的性能差异；
- 统计项目源码、汇编代码和脚本规模。

其中：

- `run_correctness.sh` 用于证明优化没有破坏协议正确性；
- `run_perf.sh` 用于证明 Montgomery 路线和 RISC-V ASM 优化带来的性能提升；
- `count_loc.sh` 用于展示项目实现规模。
