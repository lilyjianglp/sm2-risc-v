# Benchmark 说明

本文档说明面向 64 位 RISC-V 平台的 SM2 密钥交换算法高效实现项目的性能测试方法、测试版本、运行命令和主要测试结果。

本项目的 benchmark 主要用于评估 SM2 底层有限域运算、标量乘法以及 SM2 密钥交换协议层的执行效率。当前性能测试重点比较四个构建版本：

| 版本名 | 算法路线 | 实现方式 | 说明 |
|---|---|---|---|
| `baseline-c` | Prime 普通模 p 路线 | 全 C | 初版 baseline，不使用 ASM |
| `baseline-asm` | Prime 普通模 p 路线 | fp ASM + fn ASM | 初版路线加入 RISC-V 汇编优化 |
| `mont-c` | Montgomery 路线 | 全 C | Montgomery C 版本 |
| `mont-asm` | Montgomery 路线 | fp ASM + fn ASM | 当前最终高性能版本 |

这样划分可以分别观察：

1. Prime 普通模 p 路线与 Montgomery 路线的性能差异；
2. C 实现与 RISC-V ASM 实现的性能差异；
3. 最终 `mont-asm` 相比初始 `baseline-c` 的整体加速效果。

---

## 1. 测试环境

当前测试主要在 SpacemiT M1 / Muse Pi Pro RISC-V 板端完成。

| 项目 | 说明 |
|---|---|
| 平台 | SpacemiT M1 / Muse Pi Pro |
| 架构 | RISC-V 64-bit |
| 编译器 | gcc |
| 性能计数 | Linux `perf`、RISC-V `rdcycle` |
| 主要测试对象 | SM2 标量乘法、SM2 KEX、有限域运算 |

在使用 Linux `perf` 前，需要开启用户态性能计数权限：

```bash
sudo sysctl -w kernel.perf_user_access=2
```
## 2. 编译与运行命令

进入 src/ 目录：
```bash
cd src
```
编译全部版本：
```bash
make clean
make all
```
运行完整 perf 测试：
```bash
sudo sysctl -w kernel.perf_user_access=2
make perf-all
```
如果需要保存测试结果，可以使用：
```bash
make perf-all 2>&1 | tee ../bench/results/perf_all.txt
```
这里需要加 2>&1，因为 perf stat 的统计结果默认输出到 stderr。

## 3. 四个版本对应关系

当前 Makefile 中四个测试程序与文档版本名对应如下：

|文档版本名	|Makefile 目标|	说明|
|--|--|--|
baseline-c	|test_prime_c|	Prime 普通模 p 路线，全 C
baseline-asm|	test_prime_fp_asm	|Prime 普通模 p 路线，fp ASM + fn ASM
mont-c|	test_mont_c	|Montgomery 路线，全 C
mont-asm	|test_mont_asm	|Montgomery 路线，fp ASM + fn ASM

make perf-all 会依次对这四个程序执行：
```bash
perf stat -e cycles,instructions ./test_prime_c
perf stat -e cycles,instructions ./test_prime_fp_asm
perf stat -e cycles,instructions ./test_mont_c
perf stat -e cycles,instructions ./test_mont_asm
```
## 4. 整体 perf 测试结果
|版本	|Cycles|	Instructions	|IPC	|Time|	相对 baseline-c
|--|--|--|--|--|--|
|baseline-c|	19.40B|	31.92B|	1.64|	12.13s|	1.00x|
|baseline-asm|	24.61B|	35.61B	|1.45	|15.38s|	0.79x
|mont-c|	11.12B|	14.32B|	1.29	|6.95s|	1.75x|
|mont-asm	|9.55B|	12.80B|	1.34	|5.97s|	2.03x|

从整体结果看，`mont-asm `是当前最快版本。相比初始` baseline-c`，运行时间从 12.13s 降低到 5.97s，整体约获得 **2.03x 加速**。

## 5. 关键路径性能结果
指标	|`baseline-c`	|`mont-c`|	`mont-asm`	|`baseline-c → mont-asm` 提升
|--|--|--|--|--|
ScalarMul	|454,053|	305,038	2|56,810	|1.77x|
KEX Initiator Gen RA	|666,089|	422,925	|360,397|	1.85x|
KEX Responder full	|6,344,837|	4,220,395	|3,527,964	|1.80x|
KEX Initiator verify|	6,333,281	|3,786,915	|3,165,830	|2.00x|

可以看到，Montgomery 路线对 SM2 密钥交换关键路径的优化非常明显。最终 `mont-asm` 在标量乘法、KEX 发起方、KEX 响应方和验证路径上都获得了稳定提升，其中 `KEX Initiator verify `阶段接近 **2.00x 加速**。

## 6. GmSSL 函数级对比

除四版本整体 perf 测试外，项目还提供 `bench_compare/`，用于将本项目实现与 GmSSL 参考实现进行函数级性能对比。

该测试主要用于观察底层核心函数的性能差异。输出中：

```text
c     = 本项目 C 实现
asm   = 本项目 RISC-V ASM 实现
gmssl = GmSSL 参考实现
```
单位：cycles / operation，数值越低表示性能越好。

|层级	|测试项	|C|	ASM|	GmSSL|	speedup ASM/GmSSL|
|--|--|--|--|--|--|
|Fp	|fp_mont_mul|	397	|324|	1535	|4.74x|
|Fp	|fp_mont_sqr	|337|	288	|1548|	5.38x|
|Fp	|fp_from_mont|	377	|156|	1545|	9.90x|
|Fp	|fp_mont_inv	|21686	|17808	|437734	|24.58x|
|Fn|	fn_mont_mul	|952|	631|	1536	|2.43x|
|Fn	|fn_mul	|3911	|2625	|6192|	2.36x|
|ScalarMul|	point_mul_G	|289875|	252862|	646270	|2.56x|
|ScalarMul	|point_mul_P|	1610542	1|349850	|4988987|	3.70x|

这组数据说明，本项目在 Montgomery 有限域运算、标量域运算和椭圆曲线标量乘法等底层核心路径上，相比 GmSSL 参考实现取得了明显性能优势。其中 `fp_mont_inv` 的提升最突出，ASM 版本相较 GmSSL 约达到 24.58x 加速。
## 7. 结果分析
### 7.1 Montgomery 路线是主要优化来源

`mont-c` 相比 `baseline-c` 已经有明显提升：
```
baseline-c: 12.13s
mont-c:      6.95s
```
说明 Montgomery 路线本身有效降低了 SM2 运算中的模运算开销。

### 7.2 mont-asm 是当前最快版本

`mont-asm `在 Montgomery 路线基础上继续使用 RISC-V ASM 优化，最终运行时间为：
```
mont-asm: 5.97s
```
相比 `baseline-c`：
```
12.13s → 5.97s
```
整体约获得：
```
2.03x speedup
```
因此，`mont-asm` 是当前项目的最终高性能版本。

### 7.3 baseline-asm 慢于 baseline-c 的说明

测试结果中，`baseline-asm` 比 `baseline-c` 慢：
```
baseline-c:    12.13s
baseline-asm:  15.38s
```
这说明 Prime 普通模 p 路线下的手写 ASM 并没有带来有效加速。可能原因包括指令数量增加、数据依赖链较长、IPC 下降，或者该 ASM 路径没有命中当前 SM2 KEX 的主要性能瓶颈。

这个结果本身是有价值的对照：它说明本项目最终性能提升主要来自 Montgomery 路线，以及在 Montgomery 路线基础上的 ASM 优化，而不是简单地“加入 ASM 就一定更快”。

## 8. IPC 说明

IPC 是 `Instructions Per Cycle`，即每个 CPU 周期平均执行多少条指令：
```
IPC = instructions / cycles
```
IPC 高不一定代表程序更快，还需要结合总 `instructions`、总 `cycles` 和总运行时间一起看。

例如：

baseline-c IPC = 1.64
mont-asm IPC   = 1.34

虽然` baseline-c` 的 IPC 更高，但它总共执行了 31.92B 条指令；而` mont-asm` 只执行了 12.80B 条指令，所以` mont-asm `最终更快。

因此，最终性能主要看总` cycles`、总 `instructions` 和 elapsed time，而不是只看 IPC。


## 9. 注意事项
benchmark 前建议先运行正确性测试：
```bash
make correctness
```
`perf` 测试前需要开启用户态性能计数权限：
```bash
sudo sysctl -w kernel.perf_user_access=2
```
不同平台、编译器版本、CPU 频率策略都会影响 benchmark 结果。
对于波动较大的测试项，可以增加迭代次数或重复运行多次取平均值。