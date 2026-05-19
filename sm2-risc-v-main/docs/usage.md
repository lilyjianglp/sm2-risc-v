# 使用说明

本文档说明面向 64 位 RISC-V 平台的 SM2 密钥交换算法高效实现项目的常用运行方式，包括编译、正确性测试、性能测试、SM4 安全通信 demo、前端展示和一键复现脚本。

如果只是快速体验项目，建议按下面顺序运行：

```bash
cd src
make clean
make all
make correctness
make perf-all
make secure_message_demo
./secure_message_demo
```
## 1. 进入项目目录

假设项目已经位于本地：
```bash
cd /path/to/SM2-RISC-V
```
核心后端代码位于：
```
src/
```
进入后端目录：
```bash
cd src
```
## 2. 编译 SM2 benchmark 目标

项目当前主要包含四个构建版本：

|版本名|	Makefile 目标	|算法路线|	实现方式
|--|--|--|--|
baseline-c|	test_prime_c|	Prime 普通模 p 路线	|全 C
baseline-asm|	test_prime_fp_asm|	Prime 普通模 p 路线	|fp ASM + fn ASM
mont-c	|test_mont_c|	Montgomery 路线|	全 C
mont-asm|	test_mont_asm|	Montgomery 路线	|fp ASM + fn ASM

编译全部版本：
```bash
make clean
make all
```
编译完成后会生成：
```
test_prime_c
test_prime_fp_asm
test_mont_c
test_mont_asm
```
## 3. 运行正确性测试

正确性测试用于验证四种实现路径在 SM2 密钥交换协议层面输出一致。

运行命令：
```bash
make correctness
```
该命令会依次运行：
```bash
./test_correct_prime_c
./test_correct_prime_fp_asm
./test_correct_mont_c
./test_correct_mont_asm
```
测试内容包括：

- 固定 SM2 KEX 测试向量；
- Alice / Bob 双方共享密钥一致性验证；
- S1 / S2 校验值验证；
- 10 轮随机 SM2 KEX 测试；
- 四个版本固定测试向量输出一致性检查。

当前四个版本固定测试向量输出一致：
```
SharedKey = db0e9fdca7c9b8cde2edf89e6ec184ffeb3c34dbe50bc8ad72c11f7735f4108f
S1        = 53752ebb074279b5c909d88c4ae9bb611e5dc90dba3eb079ce1133d62ff489d3
S2        = 431f7b64ceaa65f8cadffcc5b6056530c2139b4f2fd26b2fb6d05e251c709b85
```
详细说明见：```docs/correctness.md```
## 4. 运行普通 benchmark

运行完整 benchmark：
```bash
make bench-all
```
如果只测试 Montgomery 路线：
```bash
make bench-mont
```
普通 benchmark 主要用于查看程序内部 ```rdcycle``` 统计的平均```cycles```。

## 5. 运行 Linux perf 性能测试

在板端先开启 perf 用户访问权限：
```bash
sudo sysctl -w kernel.perf_user_access=2
```
运行完整 perf 测试：
```bash
make perf-all
```
如果需要保存结果：
```bash
make perf-all 2>&1 | tee ../bench/results/perf_all.txt
```
这里需要加 2>&1，因为 perf stat 的统计结果默认输出到 stderr。

也可以单独测试某个版本：
```bash
make perf-prime-c
make perf-prime-asm
make perf-mont-c
make perf-mont-asm
```
当前整体性能结果摘要：

|版本|	Time|	相对 baseline-c
|--|--|--
baseline-c|	12.13s|	1.00x
baseline-asm	|15.38s|	0.79x
mont-c|	6.95s|	1.75x
mont-asm	|5.97s|	2.03x

完整性能结果见：`docs/benchmark.md`
## 6. 运行 GmSSL 对比 benchmark

如果需要和 GmSSL 做函数级性能对比，进入：
```bash
cd ../bench_compare
```
运行：
```bash
make clean
make
make run
make perf
```
注意：`bench_compare/Makefile` 依赖本地 GmSSL 路径配置。如果本地路径不同，需要修改：

- SM2_SRC
- GMSSL_INC
- GMSSL_LIB

`bench_compare` 输出中：
```
c     = 本项目 C 实现
asm   = 本项目 RISC-V ASM 实现
gmssl = GmSSL 参考实现
```
该测试主要用于比较底层函数，例如：
```
fp_mont_mul
fp_mont_sqr
fp_mont_inv
fn_mont_mul
fn_mul
point_mul_G
point_mul_P
```
## 7. 运行 SM4 安全通信 demo

SM4 demo 用于展示 SM2 密钥交换结果在安全通信中的实际用途。

回到`src/` 目录：
```bash
cd ../src
```
编译并运行：
```bash
make secure_message_demo
./secure_message_demo
```
CLI 菜单如下：
```
=========================================
 RV-GM-Secure: RISC-V GM Crypto Backend
=========================================
1. Run SM2 key exchange demo
2. Run SM4 secure message demo
3. Run correctness test
4. Run performance benchmark
0. Exit
=========================================
```
选择 `2` 后，输入明文消息，程序会自动完成：
```
输入明文
→ 运行 optimized SM2 KEX
→ 检查 KA == KB
→ KDF 派生 SM4 会话密钥
→ 生成 IV
→ SM4-CBC 加密
→ SM4-CBC 解密
→ 判断解密结果是否等于原始明文
→ 输出 PASS / FAIL
```
示例：
```
Please input message:
hello im Aria

[SM2 KEX] Shared key generated successfully.
[SM2 KEX] KA == KB: true
[SM4] Plaintext: hello im Aria
[PASS] Secure message demo success.
```
## 8. 启动前端展示页面

前端位于：`frontend/`

在 Windows 环境下可执行：
```bash
cd C:\Users\Desktop\SM2\frontend
$env:Path = "C:\Program Files\nodejs;$env:Path"
npm.cmd run dev -- --host 127.0.0.1 --port 3000
```
浏览器打开：
```
http://127.0.0.1:3000/
```
前端不是密码算法实现层，也不直接连接 RISC-V 板端程序。它主要用于展示：

- SM2 优化链路；
- SM2 KEX 结果；
- SM4 安全通信流程；
- RISC-V 后端调用链；
- 板端 benchmark 数据；
- CLI 运行指南。

前端样例数据位于：`frontend/src/demoData.js`
如果后续在板端跑出新的 CLI 结果，可以继续追加到 `demoData.js`。

## 9. 一键复现脚本

项目提供 `scripts/` 目录，用于快速复现测试流程：
```
scripts/
├── run_correctness.sh
├── run_perf.sh
└── count_loc.sh
```
### 9.1 正确性测试
```bash
bash scripts/run_correctness.sh
```
该脚本会进入 `src/` 目录，执行：
```bash
make clean
make correctness
```
### 9.2 性能测试
```bash
bash scripts/run_perf.sh
```
该脚本会进入 `src/` 目录，开启 perf 用户访问权限，并执行：
```bash
make clean
make all
make perf-all
```
### 9.3 代码行数统计
```bash
bash scripts/count_loc.sh
```
该脚本用于统计 `src/` 和 `scripts/` 目录下的 C 源码、头文件、汇编文件、脚本和 Makefile 行数。

## 10. 推荐演示流程

现场演示时，建议按下面顺序：
```bash
cd src
make clean
make all
make correctness
make perf-all
make secure_message_demo
./secure_message_demo
```
演示重点：

- `make correctness`：证明四种实现路径输出一致；
- `make perf-all`：展示 `mont-asm` 相比 `baseline-c` 的整体加速；
- `./secure_message_demo`：展示 SM2 shared key 进一步用于 SM4-CBC 安全通信；
- 前端页面：展示项目结构、优化链路和板端运行结果。
