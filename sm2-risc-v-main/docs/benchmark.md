# 性能测试说明

本项目使用 Linux `perf stat` 对四个 SM2 构建版本进行性能测试，主要统计 `cycles`、`instructions`、`IPC` 和运行时间。

测试目标是比较不同实现路径下 SM2 密钥交换的性能差异，并分析 Montgomery 路线和 RISC-V ASM 优化带来的效果。

---

## 1. 测试平台

```text
开发板：SpacemiT K1 / Muse Pi Pro RISC-V board
CPU 调优参数：spacemit-x60
测试工具：Linux perf
```

---

## 2. 测试版本

本次性能测试覆盖四个版本：

| 版本 | 算法路线 | 实现方式 | 说明 |
|---|---|---|---|
| `baseline-c` | Prime 普通模 p 路线 | fp C + fn C | 初版全 C baseline |
| `baseline-asm` | Prime 普通模 p 路线 | fp ASM + fn ASM | 初版路线加入 RISC-V ASM |
| `mont-c` | Montgomery 路线 | fp C + fn C | Montgomery 优化 C 版本 |
| `mont-asm` | Montgomery 路线 | fp ASM + fn ASM | Montgomery + RISC-V ASM 最终优化版本 |

---

## 3. 测试命令

在开发板上先开启 `perf` 用户访问权限：

```bash
sudo sysctl -w kernel.perf_user_access=2
```

然后进入 `src/` 目录执行：

```bash
make clean
make all
make perf-all
```

也可以在项目根目录下使用一键脚本：

```bash
bash scripts/run_perf.sh
```

---

## 4. perf-all 执行内容

`make perf-all` 会分别对四个版本执行：

```bash
perf stat -e cycles,instructions ./test_prime_c
perf stat -e cycles,instructions ./test_prime_fp_asm
perf stat -e cycles,instructions ./test_mont_c
perf stat -e cycles,instructions ./test_mont_asm
```

对应关系如下：

| 可执行文件 | 对应版本 |
|---|---|
| `test_prime_c` | `baseline-c` |
| `test_prime_fp_asm` | `baseline-asm` |
| `test_mont_c` | `mont-c` |
| `test_mont_asm` | `mont-asm` |

---

## 5. 性能测试结果

| 版本 | Cycles | Instructions | IPC | Time | 相对 baseline-c |
|---|---:|---:|---:|---:|---:|
| `baseline-c` | 19.40B | 31.92B | 1.64 | 12.13s | 1.00x |
| `baseline-asm` | 24.61B | 35.61B | 1.45 | 15.38s | 0.79x |
| `mont-c` | 11.12B | 14.32B | 1.29 | 6.95s | 1.75x |
| `mont-asm` | 9.55B | 12.80B | 1.34 | 5.97s | 2.03x |

---

## 6. 指标说明

### cycles

`cycles` 表示 CPU 执行程序消耗的时钟周期数。

一般来说，`cycles` 越少，程序执行开销越低。

### instructions

`instructions` 表示程序执行过程中完成的指令数量。

指令数越少，通常说明实现路径更精简。

### IPC

`IPC` 是 Instructions Per Cycle，即每个 CPU 周期平均执行多少条指令。

计算方式为：

```text
IPC = instructions / cycles
```

需要注意的是，`IPC` 高不一定代表程序更快。最终性能需要结合 `cycles`、`instructions` 和 `elapsed time` 一起判断。

例如，`baseline-c` 的 IPC 高于 `mont-asm`，但 `baseline-c` 的总指令数和总 cycles 明显更多，因此最终运行时间更长。

### Time

`Time` 表示程序整体运行时间。

本文中速度提升主要以 `baseline-c` 为基准进行比较。

---

## 7. 结果分析

### 7.1 Montgomery 路线效果明显

`mont-c` 相比 `baseline-c` 已经有明显提升：

```text
baseline-c: 12.13s
mont-c:      6.95s
```

这说明 Montgomery 路线本身有效降低了 SM2 密钥交换中的模运算开销。

从 cycles 看：

```text
baseline-c: 19.40B cycles
mont-c:     11.12B cycles
```

整体约获得：

```text
1.75x speedup
```

---

### 7.2 mont-asm 是当前最快版本

`mont-asm` 是当前性能最好的版本：

```text
mont-asm: 5.97s
```

相比 `baseline-c`：

```text
12.13s -> 5.97s
```

整体约获得：

```text
2.03x speedup
```

从 cycles 看：

```text
baseline-c: 19.40B cycles
mont-asm:    9.55B cycles
```

说明 Montgomery 表示和 RISC-V ASM 优化结合后，能够进一步降低执行开销。

---

### 7.3 baseline-asm 未取得加速

测试结果中，`baseline-asm` 比 `baseline-c` 更慢：

```text
baseline-c:    12.13s
baseline-asm:  15.38s
```

这说明 Prime 路线下的 ASM 优化没有在当前平台上形成有效加速。

可能原因包括：

- Prime 路线下的 ASM 实现没有命中主要性能瓶颈；
- 手写汇编指令调度不一定优于编译器生成代码；
- 大整数运算中的数据依赖链较长；
- 当前 SpacemiT X60 微架构更适合 Montgomery 路线下的优化方式。

因此，`baseline-asm` 作为对照版本保留，用于分析不同优化策略的效果。

---

## 8. 结论

性能测试结果表明：

- Montgomery 路线是主要性能提升来源；
- `mont-c` 相比 `baseline-c` 已经获得明显加速；
- `mont-asm` 是当前最快版本；
- `mont-asm` 相比 `baseline-c` 总运行时间从 12.13s 降低到 5.97s，整体约获得 2.03x 加速；
- Prime 路线下的 ASM 优化在当前平台上未取得加速，因此最终展示版本应以 `mont-asm` 为主。
