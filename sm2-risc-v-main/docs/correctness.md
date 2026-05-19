# 正确性验证

本文档说明面向 64 位 RISC-V 平台的 SM2 密钥交换算法高效实现项目的 SM2 KEX 正确性测试设计、测试命令、测试内容和测试结果。

本项目的目标是在不改变 SM2 密钥交换协议输出结果的前提下，对底层 Montgomery 运算、标量乘法和 RISC-V 汇编路径进行优化。因此，在进行性能 benchmark 之前，需要先验证不同实现路径在协议层输出是否一致。

---

## 1. 测试入口

正确性测试文件位于：

```text
src/test_correct.c
```
测试命令如下：
```bash
cd src
make clean
make correctness
```
`make correctness` 会自动编译并依次运行四个正确性测试程序。

## 2. 测试设计原则

本项目的正确性测试采用“分版本独立编译、独立运行”的设计。

这样做的原因是：不再把 normal 和 Montgomery 两套源码链接进同一个可执行文件，从而避免同名函数导致的 `multiple definition` 问题。每个版本单独编译，并运行当前实现路径下的 SM2 KEX 正确性测试，最后再比较四个版本在固定测试向量下的输出是否一致。

这种设计可以分别验证：

- Prime 普通模 p 路线和 Montgomery 路线是否一致；
- C 实现和 RISC-V ASM 实现是否一致；
- 底层优化是否影响 SM2 KEX 协议层输出。
## 3. 测试版本

`make correctness` 会依次运行以下四个程序：
```bash
./test_correct_prime_c
./test_correct_prime_fp_asm
./test_correct_mont_c
./test_correct_mont_asm
```
对应关系如下：

|测试程序	|版本名称|	说明|
|--|--|--|
test_correct_prime_c	|baseline-c	|普通模 p 路线，C 实现
test_correct_prime_fp_asm	|baseline-asm|	普通模 p 路线，fp 层 ASM 实现
test_correct_mont_c	|mont-c|	Montgomery 路线，C 实现
test_correct_mont_asm	|mont-asm	|Montgomery 路线，fp 层 ASM 实现

其中：

- baseline-c 用作普通 C 基准实现；
- baseline-asm 用于验证普通路线下 fp 层 ASM 优化不会破坏结果；
- mont-c 用于验证 Montgomery 路线的 C 实现正确性；
- mont-asm 用于验证 Montgomery 路线下 fp 层 ASM 优化后的协议输出正确性。
## 4. 测试内容

每个版本都会执行 SM2 密钥交换正确性测试，主要包括以下内容：

1. 固定 SM2 KEX 测试向量
2. Alice / Bob 双方共享密钥一致性验证
3. S1 / S2 校验值验证
4. 10 轮随机 SM2 KEX 测试
5. 四个版本固定测试向量输出一致性检查
### 4.1 固定 SM2 KEX 测试向量

固定测试向量用于保证每次测试输入一致，便于比较不同实现路径的输出结果。

测试会检查：

- Alice 是否能成功生成共享密钥；
- Bob 是否能成功生成共享密钥；
- 双方共享密钥是否一致；
- S1 和 S2 校验值是否正确。
### 4.2 Alice / Bob 共享密钥一致性验证

SM2 KEX 的核心目标是让通信双方在不直接传输会话密钥的情况下，计算得到一致的共享密钥。

测试会验证：
```
KA == KB
```
其中：

|符号|	含义|
|--|--|
KA	|Alice / Initiator 侧计算出的共享密钥
KB	|Bob / Responder 侧计算出的共享密钥

如果 `KA == KB`，说明双方密钥协商成功。

### 4.3 S1 / S2 校验值验证

除了共享密钥一致性外，SM2 KEX 还需要通过确认值验证双方是否确实计算出了相同的中间结果。

测试会检查：
```
S1
S2
```
其中：

|校验值|	说明|
|--|--|
S1	|一方用于验证对端密钥确认结果
S2	|另一方返回的密钥确认结果

如果 `S1 / S2` 校验通过，说明密钥交换过程中的确认值计算也保持一致。

### 4.4 随机 SM2 KEX 测试

除固定测试向量外，每个版本还会执行 10 轮随机 SM2 KEX 测试，用于验证实现不是只对固定输入正确，而是在不同随机私钥和临时密钥下仍然能够保持协议正确性。

测试结果形式为：
```
10/10 PASS
```
表示 10 轮随机测试全部通过。

### 4.5 四个版本输出一致性检查

最后，测试会检查四个版本在固定测试向量下的输出是否完全一致，包括：
```
SharedKey
S1
S2
```
这一步用于证明普通 C、普通 ASM、Montgomery C、Montgomery ASM 四种实现路径在协议层输出一致。

## 5. 测试结果

当前四个版本均通过正确性测试。

|版本|	固定 KEX 测试	|随机 KEX 测试	|结果|
|--|--|--|--|
baseline-c	|PASS	|10/10 PASS|	PASS
baseline-asm	|PASS|	10/10 PASS|	PASS
mont-c|	PASS|	10/10 PASS	|PASS
mont-asm|	PASS|	10/10 PASS|	PASS

四个版本固定测试向量输出完全一致：
```
SharedKey = db0e9fdca7c9b8cde2edf89e6ec184ffeb3c34dbe50bc8ad72c11f7735f4108f
S1        = 53752ebb074279b5c909d88c4ae9bb611e5dc90dba3eb079ce1133d62ff489d3
S2        = 431f7b64ceaa65f8cadffcc5b6056530c2139b4f2fd26b2fb6d05e251c709b85
```
## 6. 结论

根据当前测试结果：
```
baseline-c、baseline-asm、mont-c、mont-asm 四种实现路径
在 SM2 密钥交换协议层面输出一致。
```
这说明：

- 普通模 p 路线和 Montgomery 路线在 SM2 KEX 协议层输出一致；
- C 实现和 RISC-V ASM 实现输出一致；
- Montgomery 优化没有改变协议计算结果；
- RISC-V 汇编优化没有破坏协议正确性；
- 当前性能 benchmark 可以建立在正确性通过的基础上进行比较。

因此，可以认为当前优化版本在 SM2 密钥交换协议层面保持正确。

## 7. 与 benchmark 的关系

正确性测试应当先于性能 benchmark 执行。

推荐流程如下：
```bash
cd src
make clean
make correctness
make bench-all
make perf-all
```
只有当 `make correctness` 通过后，后续 benchmark 结果才具有比较意义。

## 8. 后续补充方向

后续可以继续完善以下正确性测试内容：

- 增加更多轮随机 SM2 KEX 测试；
- 增加边界私钥和边界临时密钥测试；
- 增加异常输入测试，例如非法点、无穷远点、零标量等；
- 增加 SM4 demo 端到端正确性测试；
- 将 correctness 测试加入 CI 流程。