# 安装说明

本文档说明 SM2 RISC-V 项目的基础运行环境、依赖安装、核心代码编译方法，以及可选的 GmSSL 对照测试环境配置方法。

SM2 RISC-V 当前主要面向 SpacemiT Muse Pi Pro / SpacemiT M1(X60) RISC-V 平台进行测试和优化。项目核心代码位于 `src/` 目录，主要依赖 `gcc`、`make`、`binutils` 和 Linux `perf`。如果需要运行与 GmSSL 的对照 benchmark，还需要提前准备 GmSSL 并配置对应路径。

---

## 1. 环境要求

推荐环境如下：

| 项目 | 说明 |
|---|---|
| 硬件平台 | SpacemiT Muse Pi Pro / SpacemiT M1(X60) |
| 架构 | RISC-V 64-bit |
| 操作系统 | Linux |
| 编译器 | gcc |
| 构建工具 | make |
| 二进制工具 | objcopy，来自 binutils |
| 性能分析工具 | Linux perf |
| 对照实现 | GmSSL，可选 |

当前 `src/Makefile` 中默认使用 SpacemiT X60 平台相关编译参数：

```makefile
CPU_TUNE   := spacemit-x60
CPU_CFLAGS := -mcpu=spacemit-x60
```
因此，建议优先在 SpacemiT Muse Pi Pro / X60 实机环境中编译和测试。

如果当前 gcc 不支持 -mcpu=spacemit-x60，编译时可能会报错。此时需要更换支持该平台的工具链，或者根据本机工具链支持情况修改 src/Makefile 中的 CPU_CFLAGS。

## 2. 安装基础依赖

在 Linux 环境中安装基础构建工具：
```bash
sudo apt update
sudo apt install -y build-essential make gcc binutils git
```
各依赖作用如下：

|依赖	|作用|
|--|--|
gcc	|编译 C 源码和 RISC-V 汇编文件
make|	执行项目中的 Makefile
binutils|	提供 objcopy 等工具，bench_compare 中会用到
git	|获取项目代码或第三方依赖
## 3. 准备 perf

项目使用 Linux perf 统计性能事件，例如：
- cycles
- instructions
  
可以先检查当前系统是否已经安装 perf：
```bash
perf --version
```
如果没有安装，需要根据当前 Linux 发行版和内核版本安装对应的 perf 工具。

在 Ubuntu / Debian 系统中，可以尝试：
```bash
sudo apt install -y linux-perf
```
如果该包不可用，可以根据系统提示安装对应版本的 linux-tools 或其他 perf 包。

## 4. 开启 perf 权限

在运行 perf stat 前，可能需要打开用户态性能计数权限。

优先尝试：
```bash

sudo sysctl -w kernel.perf_user_access=2
```
如果当前内核不支持 kernel.perf_user_access，可以尝试：
```bash
sudo sysctl -w kernel.perf_event_paranoid=0
```
或者：
```bash
sudo sysctl -w kernel.perf_event_paranoid=1
```
不同 Linux 内核和发行版的权限控制项可能不同，具体以当前系统支持情况为准。

## 5. 获取项目代码

如果项目已经在本地，进入项目根目录：

```bash
cd /path/to/SM2-RISC-V
```
如果从 Git 仓库获取：
```
git clone <your-repo-url>
cd SM2-RISC-V
```
项目主要目录如下：
```
SM2-RISC-V/
├── src/
├── src/asm/
├── include/
├── sm4_secure_demo/
├── bench_compare/
├── frontend/
├── scripts/
├── docs/
└── README.md
```
其中：

|目录|	说明|
|--|--|
src/	|SM2 核心实现、RISC-V 汇编优化、benchmark 和 CLI 后端入口
src/asm/	|RISC-V 汇编优化文件
include/|	SM4 通信模块头文件
sm4_secure_demo/	|SM4 安全通信 demo 的 CLI 展示入口和说明文档
bench_compare/	|与 GmSSL 进行统一对照 benchmark
frontend/|	React 前端展示页面
scripts/	|一键正确性测试、性能测试和代码行数统计脚本
docs/	|项目文档
## 6. 编译核心代码

核心编译入口位于 src/ 目录。
```bash
cd src
make clean
make all
```
根据当前 src/Makefile，make all 会生成以下四个目标程序：

|程序	|说明|
|--|--|
test_prime_c|	Prime / 普通模 p 路线，C 实现
test_prime_fp_asm|	Prime / 普通模 p 路线，fp 层 ASM 实现
test_mont_c|	Montgomery 路线，C 实现
test_mont_asm|	Montgomery 路线，fp 层 ASM 实现
## 7. 运行正确性测试

在进入性能测试前，建议先运行正确性测试：

```bash
make correctness
```
该命令会依次运行四个版本的 SM2 KEX 正确性测试：
```
test_correct_prime_c
test_correct_prime_fp_asm
test_correct_mont_c
test_correct_mont_asm
```
用于验证 `baseline-c`、`baseline-asm`、`mont-c`、`mont-asm` 四种实现路径在固定测试向量和随机 KEX 测试下输出一致。
## 8. 运行基础 benchmark

在 src/ 目录下运行完整 benchmark：
```bash
make bench-all
```
如果只测试 Prime 路线：
```bash
make bench-prime
```
如果只测试 Montgomery 路线：
```bash
make bench-mont
```
## 9. 使用 perf 进行性能测试

在 src/ 目录下运行完整 perf 测试：
```bash
make perf-all
```
也可以单独运行某一类测试：
```bash
make perf-prime-c
make perf-prime-asm
make perf-mont-c
make perf-mont-asm
```
这些目标分别对应：

|命令	|说明|
|--|--|
make perf-prime-c|	Prime / 普通模 p 路线，C 实现
make perf-prime-asm	|Prime / 普通模 p 路线，fp 层 ASM 实现
make perf-mont-c	|Montgomery 路线，C 实现
make perf-mont-asm|	Montgomery 路线，fp 层 ASM 实现
make perf-all	|依次运行上述所有 perf 测试
## 10. 运行专项 perf workload

src/Makefile 还提供了基于 test_perf.c 的专项 perf workload，用于分析 SM2 KEX 相关关键路径。

可用命令包括：
```bash
make perf-test-perf-responder
make perf-test-perf-initiator
make perf-test-perf-scalar
make perf-test-perf-genra
```
对应含义如下：

|命令	|说明|
|--|--|
make perf-test-perf-responder	|分析 SM2 KEX responder 路径
make perf-test-perf-initiator	|分析 SM2 KEX initiator 路径
make perf-test-perf-scalar	|分析固定基点标量乘法路径
make perf-test-perf-genra|	分析临时公钥生成路径

如果只想直接运行 test_perf：
```bash
make run-test-perf
```
## 11. 编译并运行 SM4 安全通信 demo

如果需要运行 SM4 安全通信演示，可以在 `src/` 目录下执行：

```bash
make secure_message_demo
./secure_message_demo
```
该 demo 会通过 CLI 菜单运行 SM2 KEX、KDF、SM4-CBC 加密和解密流程，用于展示 `SM2 shared key` 在上层安全通信中的实际用途。

如果当前 Makefile 中尚未包含 `secure_message_demo` 目标，需要先确认 `src/`Makefile 是否已经合并了 SM4 demo 的编译规则。

## 12. GmSSL 对照测试，可选

如果只运行 src/ 下的 benchmark，可以暂时不配置 GmSSL。

如果需要运行 bench_compare/ 中的 GmSSL 对照 benchmark，则需要提前准备 GmSSL，并配置 bench_compare/Makefile 中的路径。

当前 bench_compare/Makefile 中涉及三个路径变量：
```bash
SM2_SRC   := /path/to/RV-GM-Secure/src
GMSSL_INC := /path/to/GmSSL/include
GMSSL_LIB := /path/to/GmSSL/build/bin
```
含义如下：

|变量	|说明|
|--|--|
SM2_SRC	|本项目 src/ 目录路径
GMSSL_INC	GmSSL| 头文件目录
GMSSL_LIB	GmSSL| 编译后库文件所在目录

例如，如果项目和 GmSSL 位于 /home/sm2/ 下，可以配置为：
```
SM2_SRC   := /home/sm2/sm2-risc-v-main/sm2-risc-v-main/src
GMSSL_INC := /home/sm2/GmSSL/include
GMSSL_LIB := /home/sm2/GmSSL/build/bin
```
## 13. 查看 bench_compare 当前配置

进入 bench_compare/ 后，可以先查看当前路径配置：
```bash
cd bench_compare
make show
```
该命令会输出：
```
SM2_SRC=...
GMSSL_INC=...
GMSSL_LIB=...
CFLAGS_COMMON=...
```
如果路径不符合本机环境，需要先修改 bench_compare/Makefile。

## 14. 编译并运行 GmSSL 对照 benchmark

进入 bench_compare/：
```bash
cd bench_compare
```
清理旧文件：
```bash
make clean
```
编译统一 benchmark 程序：
```bash
make
```
运行 benchmark：
```bash
make run
```
使用 perf 统计：
```bash
make perf
```
make run 和 make perf 会通过 LD_LIBRARY_PATH=$(GMSSL_LIB) 指定 GmSSL 动态库路径。

## 15. bench_compare 的工作方式

bench_compare/Makefile 会分别编译本项目的 C 版本对象和 ASM 版本对象，然后通过 objcopy --prefix-symbols 对符号进行前缀化处理。

其中：

|前缀	|含义|
|--|--|
c_	|本项目 C 实现
asm_	|本项目 ASM 实现

例如：
```bash
$(OBJCOPY) --prefix-symbols=c_   ...
$(OBJCOPY) --prefix-symbols=asm_ ...
```
这样做的目的是避免 C 版本和 ASM 版本函数名冲突，从而能够在同一个 bench_all 程序中同时链接：

- 本项目 C 实现；
- 本项目 ASM 实现；
- GmSSL 参考实现。

最终生成的程序为：
```
bench_all
```
## 16. 动态库路径

如果运行 GmSSL 对照测试时出现动态库找不到的问题，例如：
```
error while loading shared libraries: libgmssl.so
```
可以手动设置：
```bash
export LD_LIBRARY_PATH=/path/to/GmSSL/build/bin:$LD_LIBRARY_PATH
```
或者临时运行：
```bash
env LD_LIBRARY_PATH=/path/to/GmSSL/build/bin ./bench_all
```
如果使用 make run 或 make perf，当前 bench_compare/Makefile 已经通过 env LD_LIBRARY_PATH=$(GMSSL_LIB) 进行设置。

## 17. 安装验证
### 17.1 验证核心编译
```bash
cd src
make clean
make all
```
如果生成以下文件，说明核心编译环境可用：
```
test_prime_c
test_prime_fp_asm
test_mont_c
test_mont_asm
```
### 17.2 验证基础 benchmark
```bash
make bench-all
```
如果能够输出 Prime / Montgomery、C / ASM 对比结果，说明基础 benchmark 可用。

### 17.3 验证 perf
```bash
make perf-all
```
如果能够输出 cycles 和 instructions，说明 perf 配置可用。

### 17.4 验证 GmSSL 对照 benchmark

```bash
cd ../bench_compare
make clean
make
make run
```
如果能够输出 c、asm、gmssl 三组结果，说明 GmSSL 对照 benchmark 环境可用。

## 18. 常见问题
### 18.1 gcc 不支持 -mcpu=spacemit-x60

如果出现：
```
gcc: error: unrecognized argument in option '-mcpu=spacemit-x60'
```
说明当前 gcc 不支持该 CPU 选项。可以考虑：

在 SpacemiT Muse Pi Pro / X60 实机环境中编译；
更换支持 SpacemiT X60 的 gcc；
根据本机工具链支持情况修改 src/Makefile 中的 CPU_CFLAGS。

例如，可临时改成当前工具链支持的 RISC-V 参数，但具体取值需要根据本机 gcc 支持情况确定。

### 18.2 perf 没有权限

如果运行 perf 时出现权限问题，可以尝试：
```bash
sudo sysctl -w kernel.perf_user_access=2
```
或：
```bash
sudo sysctl -w kernel.perf_event_paranoid=0
```
### 18.3 bench_compare 找不到本项目源码

如果 bench_compare 编译时报找不到 fp_montgomery.c、fn_Montgomery.c 等文件，需要检查：
```
SM2_SRC
```
是否正确指向本项目的 src/ 目录。

### 18.4 bench_compare 找不到 GmSSL 头文件

如果报错找不到：
```
gmssl/sm2.h
gmssl/sm2_z256.h
```
需要检查：
```
GMSSL_INC
```
是否正确指向 GmSSL 的 include 目录。

### 18.5 bench_compare 找不到 GmSSL 动态库

如果运行时报动态库加载失败，需要检查：
```
GMSSL_LIB
```
并设置：
```bash
export LD_LIBRARY_PATH=/path/to/GmSSL/build/bin:$LD_LIBRARY_PATH
```
### 18.6 objcopy 前缀化后链接失败

bench_compare/Makefile 中已经使用：
```
-fno-stack-protector
```
这是为了避免 `objcopy --prefix-symbols` 对` __stack_chk_* `等符号前缀化后导致链接失败。

如果手动修改编译参数后出现相关链接问题，可以检查是否误删了`-fno-stack-protector。`

## 18. 下一步

完成安装和环境配置后，可以继续阅读：

`docs/build.md`：编译说明；
`docs/correctness.md`：正确性验证；
`docs/benchmark.md`：性能测试与结果分析；
`docs/platform.md`：测试平台说明。