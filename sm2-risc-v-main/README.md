# 面向 64 位 RISC-V 平台的 SM2 密钥交换算法高效实现

## 1. 项目简介

该项目是一个面向 RISC-V 平台的国密算法优化与安全通信演示项目。项目核心工作围绕 SM2 密钥交换中的有限域 Montgomery 运算、标量乘法、Pornin 模逆和 RISC-V 汇编优化展开，并在此基础上实现了 SM4 安全通信演示模块和前端展示页面。

项目整体链路如下：

```text
SM2 密钥交换
→ Montgomery 有限域运算后端
→ Pornin 模逆与 RISC-V 汇编优化
→ 生成双方一致 shared key
→ KDF 派生 SM4 会话密钥
→ SM4-CBC 加密与解密消息
→ CLI / 前端展示运行结果
```
项目提供：

- C baseline 实现；
- RISC-V 汇编优化实现；
- GmSSL 对照验证；
- rdcycle 微基准测试；
- Linux perf 性能剖析；
- SM4-CBC 安全通信演示；
- CLI 菜单式交互；
- React 前端展示页面。

## 2. 项目背景

SM2 是国密体系中的重要公钥密码算法，广泛用于数字签名、密钥交换和身份认证等场景。在 SM2 密钥交换过程中，性能开销主要来自椭圆曲线标量乘法，而标量乘法又依赖大量底层有限域运算。

本项目关注以下问题：

- SM2 依赖大量有限域加减乘、平方和求逆运算；
- Montgomery 乘法、平方、模逆是性能关键路径；
- RISC-V 平台缺少传统硬件进位标志，需要针对指令集和流水线特点进行优化；
- 底层优化需要能够传导到 SM2 密钥交换协议层；
- 密钥交换结果需要通过上层安全通信 demo 展示其实际用途。

因此，本项目不仅实现底层性能优化，还进一步通过 SM4-CBC 安全通信 demo 展示 SM2 shared key 在会话加密中的应用。
## 3. 核心特性

- 支持 SM2 密钥交换流程；
- 支持普通模 p 路线与 Montgomery 路线；
- 支持 C 实现与 RISC-V ASM 实现对比；
- 支持 fp / fn 层底层运算 benchmark；
- 支持标量乘法与 KEX 级别 benchmark；
- 支持四版本 SM2 KEX 正确性测试；
- 支持与 GmSSL 进行性能对比；
- 支持 Linux perf 进行专项性能剖析。
- 支持基于 SM2 shared key 的 SM4-CBC 安全通信演示；
- 支持 CLI 菜单式交互；
- 支持前端展示 SM2 优化链路、SM4 安全信道和板端实测数据。

## 4. 测试平台
当前主要测试平台如下：
- SpacemiT Muse Pi Pro / X60；
- RISC-V 64-bit；
- gcc；
- Linux perf；
- RISC-V rdcycle。
  


## 5. 项目结构

项目当前采用“前端展示层 + RISC-V 后端实现层 + 文档说明”的结构。根目录下包含总 README、前端工程和 SM2 RISC-V 后端工程。

简化目录结构如下：

```text
.
├── README.md
├── docs/
│   ├── install.md
│   ├── usage.md
│   ├── benchmark.md
│   ├── correctness.md
│   ├── sm2_kex.md
│   └── sm4_demo.md
│
├── frontend/
│   ├── index.html
│   ├── package.json
│   ├── package-lock.json
│   ├── vite.config.js
│   ├── README.md
│   └── src/
│       ├── main.jsx
│       ├── styles.css
│       └── demoData.js
│
└── sm2-risc-v-main/
    ├── README.md
    ├── bench_compare/
    │   ├── bench_all.c
    │   ├── bench_compare.c
    │   ├── correctness_compare.c
    │   └── Makefile
    │
    └── src/
        ├── Makefile
        ├── fp.c / fp_montgomery.c
        ├── fn.c / fn_Montgomery.c
        ├── sm2_curve.c / sm2_curve_mont.c
        ├── sm2_scalar.c / sm2_scalar_mont.c
        ├── sm2_kex.c / sm2_kex_mont.c
        ├── sm4_channel.cpp
        ├── sm4_channel.h
        ├── test.c
        ├── test_correct.c
        ├── test_perf.c
        ├── asm/
        │   ├── fp_mont_mul.S
        │   ├── fp_mont_sqr.S
        │   ├── fp_mont_to_from.S
        │   ├── fn_mont.S
        │   └── ...
        ├── demo/
        │   ├── README.md
        │   ├── secure_message_demo.cpp
        │   └── image/
        └── test/
            ├── test_add.c
            ├── test_mul.c
            ├── test_sqr.c
            └── test_fp_mont_correctness.c
```
各目录作用如下：

目录 / 文件|	作用
|--|--|
README.md	|项目总说明，介绍整体功能、快速开始、测试结果和贡献分工
docs/	|项目文档目录，存放安装说明、使用说明、benchmark、正确性验证和模块设计文档
frontend/|	React 前端展示页面，用于展示 SM2 优化链路、SM4 安全通信流程和板端实测数据
frontend/src/main.jsx	|前端页面内容和交互逻辑
frontend/src/styles.css	|前端页面样式
frontend/src/demoData.js|	前端展示用的板端 CLI 样例数据
sm2-risc-v-main/	|RISC-V 后端核心工程
sm2-risc-v-main/bench_compare/	|与 GmSSL 进行性能对比的 benchmark 工程
sm2-risc-v-main/src/|	SM2 核心实现、RISC-V 汇编优化、SM4 demo 后端和测试入口
sm2-risc-v-main/src/asm/|	RISC-V 汇编优化文件
sm2-risc-v-main/src/demo/|	CLI 交互演示入口，包含 secure_message_demo.cpp
sm2-risc-v-main/src/test/	|fp 层基础正确性测试代码

## 6. 快速开始

本项目的主要编译入口位于 `src/Makefile`。当前项目将 SM2 实现划分为四个构建版本：`baseline-c`、`baseline-asm`、`mont-c`、`mont-asm`，用于比较 Prime 普通模 p 路线、Montgomery 路线以及 C / ASM 实现的性能差异。

| 版本名 | Makefile 目标 | 算法路线 | 实现方式 |
|---|---|---|---|
| `baseline-c` | `test_prime_c` | Prime 普通模 p 路线 | 全 C |
| `baseline-asm` | `test_prime_fp_asm` | Prime 普通模 p 路线 | fp ASM + fn ASM |
| `mont-c` | `test_mont_c` | Montgomery 路线 | 全 C |
| `mont-asm` | `test_mont_asm` | Montgomery 路线 | fp ASM + fn ASM |

常用命令如下：

```bash
cd src

# 编译主要测试目标
make clean
make all

# 正确性测试
make correctness

# 普通 benchmark
make bench-all

# perf 性能测试
sudo sysctl -w kernel.perf_user_access=2
make perf-all

# SM4 安全通信 demo
make secure_message_demo
./secure_message_demo
```
GmSSL 对比测试：
```bash
cd ../bench_compare
make clean
make
make run
make perf
```
前端展示页面：
```bash
cd <project-root>\frontend
$env:Path = "C:\Program Files\nodejs;$env:Path"
npm.cmd run dev -- --host 127.0.0.1 --port 3000
```
浏览器打开：
```
http://127.0.0.1:3000/
```
更完整的运行方式、CLI 菜单和演示流程见：`docs/usage.md`

## 7. 核心模块说明
### 7.1 fp 素域运算模块

fp 模块实现 SM2 曲线素域上的基础运算，包括模加、模减、模乘、平方、约减、求逆以及 Montgomery 表示转换，是椭圆曲线点运算的底层支撑。

### 7.2 fn 标量域运算模块

fn 模块实现 SM2 曲线阶 n 上的标量运算，主要用于私钥、随机数和标量乘法相关计算。

### 7.3 Montgomery 运算模块

Montgomery 模块用于优化大整数模乘和模约减，将频繁出现的模乘运算转换为更适合底层实现优化的形式。项目重点测试和优化：

- fp_mont_mul
- fp_mont_sqr
- fp_mont_inv
- fp_to_mont
- fp_from_mont
  
### 7.4 RISC-V 汇编优化模块

src/asm/ 目录包含 RISC-V 平台下的关键运算汇编实现，主要覆盖 Montgomery 模加、模减、模乘、模平方、Montgomery 表示转换和标量域 Montgomery 乘法等热点函数。

代表性文件包括：
```
asm/fp_mont_mul.S
asm/fp_mont_sqr.S
asm/fp_mont_to_from.S
asm/fn_mont.S
asm/pornin_full_inv.S
```
### 7.5 SM2 密钥交换模块

sm2_kex.c 和 sm2_kex_mont.c 分别实现普通路线和 Montgomery 路线下的 SM2 密钥交换流程，用于评估底层优化对协议层性能的影响。

详细说明见：`docs/sm2_kex.md`
### 7.6 SM4 安全通信演示模块

SM4 demo 用于展示 SM2 密钥交换结果在安全通信中的实际用途。

其核心流程为：
```text
SM2 KEX 输出 shared key
        ↓
SM3/KDF 派生 16 字节 SM4 session key
        ↓
SM4-CBC 加密用户输入消息
        ↓
SM4-CBC 解密密文
        ↓
验证解密结果是否等于原始明文
```
该模块不直接调用底层 RISC-V 汇编代码，也不重新实现 SM2。它只调用上层 SM2 密钥交换接口，保持了清晰的分层结构。

详细说明见：`docs/sm4_demo.md`
### 7.7 前端展示模块

前端位于 `frontend/` 目录，用于展示项目结构、SM2 优化链路、SM4 安全通信流程和 RISC-V 板端实测性能数据。需要说明的是，前端不是密码算法实现层，也不直接连接 RISC-V 板端程序；它主要作为项目展示页面，将 CLI 真实运行结果和 benchmark 数据进行可视化，方便答辩和演示。

前端页面主要包含 6 个部分：

| 页面模块 | 说明 |
|---|---|
| SM2 优化概览 | 展示项目核心工作，包括 Montgomery 有限域后端、Pornin 模逆优化、RISC-V 汇编优化和 SM2 KEX 性能提升 |
| SM2 密钥交换 | 展示 SM2 KEX 的运行结果，包括 shared key 摘要、`KA == KB` 一致性验证和板端实测状态 |
| SM4 安全信道 | 展示基于 SM2 shared key 的 SM4-CBC 加密通信流程，包括密钥派生、加密、解密和 PASS/FAIL 验证 |
| RISC-V 后端链路 | 展示从 CLI / 前端展示层到 SM2 KEX、标量乘法、Montgomery 后端、Pornin 模逆和 SM4-CBC 的整体调用链 |
| 性能基准评测 | 展示 RISC-V 板端实测性能数据，如 Montgomery 模乘、模平方、Pornin 模逆、SM2 KEX 发起方和响应方 cycles |
| 运行指南 | 展示 CLI 编译、运行、correctness test 和 benchmark 命令，便于现场复现 |

前端主要文件说明如下：

| 文件 | 作用 |
|---|---|
| `frontend/src/main.jsx` | 页面内容和交互逻辑 |
| `frontend/src/styles.css` | 页面样式 |
| `frontend/src/demoData.js` | 板端 CLI 样例数据 |
| `frontend/package.json` | 前端依赖和启动脚本 |
| `frontend/index.html` | React 挂载入口 |
| `frontend/vite.config.js` | Vite 开发服务器配置 |

前端样例数据保存在 `frontend/src/demoData.js`，当前样例来自 RISC-V 板端 CLI 实际运行输出，例如：

```text
hello im Aria
risc-v sm2 kex ready
secure channel message
SM4 demo on Muse Pi Pro
```
## 8. 正确性验证

本项目提供独立的 SM2 KEX 正确性测试入口，用于验证不同实现路径下的密钥交换结果是否一致。正确性测试文件位于：

```text
src/test_correct.c
```
测试命令如下：
```bash
cd src
make clean
make correctness
```
当前正确性测试不再将 normal 和 Montgomery 两套源码链接到同一个可执行文件中，而是分别编译并运行四个版本，避免函数符号冲突和 multiple definition 问题。

make correctness 会依次运行：
```bash
./test_correct_prime_c
./test_correct_prime_fp_asm
./test_correct_mont_c
./test_correct_mont_asm
```
四个版本分别对应：

|版本	|说明|
|--|--|
baseline-c	|普通模 p 路线，C 实现
baseline-asm	|普通模 p 路线，fp 层 ASM 实现
mont-c	|Montgomery 路线，C 实现
mont-asm|	Montgomery 路线，fp 层 ASM 实现

测试内容包括：

- 固定 SM2 KEX 测试向量；
- Alice / Bob 双方共享密钥一致性验证；
- S1 / S2 校验值验证；
- 10 轮随机 SM2 KEX 测试；
- 四个版本固定测试向量输出一致性检查。

当前四个版本均已通过正确性测试，固定测试向量输出完全一致：
```
SharedKey = db0e9fdca7c9b8cde2edf89e6ec184ffeb3c34dbe50bc8ad72c11f7735f4108f
S1        = 53752ebb074279b5c909d88c4ae9bb611e5dc90dba3eb079ce1133d62ff489d3
S2        = 431f7b64ceaa65f8cadffcc5b6056530c2139b4f2fd26b2fb6d05e251c709b85
```
这说明 baseline-c、baseline-asm、mont-c、mont-asm 四种实现路径在 SM2 密钥交换协议层面输出一致，Montgomery 优化和 RISC-V ASM 优化没有破坏协议正确性。

详细测试说明见：```docs/correctness.md```

## 9. Benchmark 结果摘要
本项目在 SpacemiT M1 / Muse Pi Pro RISC-V 板端进行性能测试。当前性能测试主要比较四个构建版本：

| 版本名 | 算法路线 | 实现方式 | 说明 |
|---|---|---|---|
| `baseline-c` | Prime 普通模 p 路线 | 全 C | 初版 baseline，不使用任何 ASM |
| `baseline-asm` | Prime 普通模 p 路线 | fp ASM + fn ASM | 初版路线加 RISC-V 汇编优化 |
| `mont-c` | Montgomery 路线 | 全 C | 优化后的 Montgomery C 版本 |
| `mont-asm` | Montgomery 路线 | fp ASM + fn ASM | 当前最终高性能版本 |


测试命令如下：

```bash
cd src
sudo sysctl -w kernel.perf_user_access=2
make clean
make perf-all
```

### 9.1 四版本整体 perf 结果
|版本|	Cycles	|Instructions|	IPC	|Time|	相对 baseline-c|
|--|--|--|--|--|--|
baseline-c	|19.40B	|31.92B	|1.64|	12.13s	|1.00x
baseline-asm	|24.61B	|35.61B	|1.45	|15.38s	|0.79x
mont-c	|11.12B	|14.32B	|1.29	|6.95s|	1.75x
mont-asm	|9.55B	|12.80B	|1.34|	5.97s	|2.03x

从整体结果看，mont-asm 是当前最快版本。相比初始 baseline-c，运行时间从 12.13s 降低到 5.97s，整体约获得 2.03x 加速。

### 9.2 关键路径性能对比
|指标|	baseline-c	|mont-c	|mont-asm	|baseline-c → mont-asm 提升|
|--|--|--|--|--|
ScalarMul	|454,053	|305,038|	256,810|	1.77x|
KEX Initiator Gen RA|	666,089|	422,925	|360,397|	1.85x|
KEX Responder full	|6,344,837|	4,220,395|	3,527,964|	1.80x
KEX Initiator verify|	6,333,281|	3,786,915	|3,165,830	|2.00x

结果表明，Montgomery 路线是主要优化来源，进一步结合 RISC-V ASM 后，SM2 标量乘法和 KEX 协议层关键路径均获得明显提升。其中 KEX Initiator verify 阶段接近 2.00x 加速。

完整 benchmark 方法、IPC 解释和更详细的结果分析见：`docs/benchmark.md`
## 10. SM4 安全通信演示

本项目已经实现 SM4 安全通信 demo，用于展示优化后的 SM2 密钥交换结果如何用于上层安全通信。

示例运行结果：
```
Please input message:
hello im Aria

[SM2 KEX] Running optimized SM2 key exchange...
[SM2 KEX] Shared key generated successfully.
[SM2 KEX] KA == KB: true

[KDF] Deriving SM4 session key from SM2 shared key...
[SM4] Encrypting message...
[SM4] Ciphertext: ef7a092e929e4f0256abb22a2ac5a1cc

[SM4] Decrypting message...
[SM4] Plaintext: hello im Aria

[PASS] Secure message demo success.
```
该结果说明：

- SM2 KEX 成功生成双方一致 shared key；
- KDF 成功派生 SM4 会话密钥；
- SM4-CBC 成功完成加密和解密；
- 解密结果与原始输入一致。
  
## 11. 一键复现脚本

项目提供 `scripts/` 目录，用于快速复现正确性测试、性能测试和代码行数统计。

```text
scripts/
├── run_correctness.sh
├── run_perf.sh
└── count_loc.sh
```
常用命令：
```bash
# 运行 SM2 KEX 正确性测试
bash scripts/run_correctness.sh

# 运行四版本 perf 性能测试
bash scripts/run_perf.sh

# 统计源码、汇编、脚本和 Makefile 行数
bash scripts/count_loc.sh
```
详细使用说明见：docs/usage.md
## 12. 贡献者分工

| 成员 | 主要职责 | 贡献内容 |
| --- | --- | --- |
| 王奕文 | 项目集成与复现流程 | 负责项目结构整理、编译流程、测试流程、benchmark 数据整理和运行脚本 |
| 姜良萍 | SM4 通信与交互演示 | 负责 SM4 通信 demo、CLI 菜单、示例程序和交互演示功能 |
| 崔钰婷 | README 与开源文档 | 负责 README、安装说明、使用说明、设计文档、benchmark 文档和开源规范建设 |