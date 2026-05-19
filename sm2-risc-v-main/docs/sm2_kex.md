
# SM2 密钥交换说明

本文档说明面向 64 位 RISC-V 平台的 SM2 密钥交换算法高效实现项目中 SM2 密钥交换模块的功能、协议流程、代码结构、普通路线与 Montgomery 路线的区别，以及相关正确性验证和性能测试入口。

SM2 密钥交换用于通信双方在不直接传输会话密钥的情况下，基于双方的静态密钥、临时密钥和身份标识协商出一致的共享密钥。该共享密钥后续可以用于 SM4 等对称加密算法，从而实现安全通信。

在本项目中，SM2 KEX 不只是一个独立协议模块，也承担了验证底层优化效果的作用。项目通过 SM2 KEX 观察有限域运算、Montgomery 运算、标量乘法和 RISC-V 汇编优化是否能够真正传导到协议层。

---

## 1. 模块位置

SM2 密钥交换相关代码主要位于 `src/` 目录：

```text
src/
├── sm2_kex.h
├── sm2_kex.c
├── sm2_kex_mont.c
├── sm2_curve.h
├── sm2_curve.c
├── sm2_curve_mont.c
├── sm2_scalar.h
├── sm2_scalar.c
├── sm2_scalar_mont.c
├── fp.c
├── fp_montgomery.c
├── fn.c
└── fn_Montgomery.c
```

其中：

| 文件 | 说明 |
|---|---|
| `sm2_kex.h` | SM2 KEX 普通版和 Montgomery 版接口声明 |
| `sm2_kex.c` | 普通域下的 SM2 KEX 实现 |
| `sm2_kex_mont.c` | Montgomery 域下的 SM2 KEX 实现 |
| `sm2_curve.c` | 普通域椭圆曲线点运算 |
| `sm2_curve_mont.c` | Montgomery 域椭圆曲线点运算 |
| `sm2_scalar.c` | 普通域标量乘法 |
| `sm2_scalar_mont.c` | Montgomery 域标量乘法 |
| `fp.c` | 普通素域运算 |
| `fp_montgomery.c` | Montgomery 素域运算 |
| `fn_Montgomery.c` | 阶 n 上的 Montgomery 标量运算 |

---

## 2. 对外接口

`sm2_kex.h` 中同时提供普通版和 Montgomery 版接口。

### 2.1 普通版接口

```c
int sm2_kex_initiator_gen_RA(sm2_affine_t *RA,
                             const uint8_t rA[32]);

int sm2_kex_initiator_compute_key(uint8_t *K, size_t klen,
                                  uint8_t S1[32],
                                  uint8_t SA[32],
                                  const uint8_t *peer_S2,
                                  const uint8_t dA[32],
                                  const sm2_affine_t *PA,
                                  const uint8_t *idA, size_t idA_len,
                                  const uint8_t rA[32],
                                  const sm2_affine_t *RA,
                                  const sm2_affine_t *PB,
                                  const uint8_t *idB, size_t idB_len,
                                  const sm2_affine_t *RB);

int sm2_kex_responder_compute_key(sm2_affine_t *RB,
                                  uint8_t *K, size_t klen,
                                  uint8_t S2[32],
                                  uint8_t SB[32],
                                  const uint8_t peer_S1[32],
                                  const uint8_t dB[32],
                                  const sm2_affine_t *PB,
                                  const uint8_t *idB, size_t idB_len,
                                  const uint8_t rB[32],
                                  const sm2_affine_t *RA,
                                  const sm2_affine_t *PA,
                                  const uint8_t *idA, size_t idA_len);
```

### 2.2 Montgomery 版接口

```c
int sm2_kex_initiator_gen_RA_mont(sm2_affine_t *RA,
                                  const uint8_t rA[32]);

int sm2_kex_initiator_compute_key_mont(uint8_t *K, size_t klen,
                                       uint8_t S1[32],
                                       uint8_t SA[32],
                                       const uint8_t *peer_S2,
                                       const uint8_t dA[32],
                                       const sm2_affine_t *PA,
                                       const uint8_t *idA, size_t idA_len,
                                       const uint8_t rA[32],
                                       const sm2_affine_t *RA,
                                       const sm2_affine_t *PB,
                                       const uint8_t *idB, size_t idB_len,
                                       const sm2_affine_t *RB);

int sm2_kex_responder_compute_key_mont(sm2_affine_t *RB,
                                       uint8_t *K, size_t klen,
                                       uint8_t S2[32],
                                       uint8_t SB[32],
                                       const uint8_t peer_S1[32],
                                       const uint8_t dB[32],
                                       const sm2_affine_t *PB,
                                       const uint8_t *idB, size_t idB_len,
                                       const uint8_t rB[32],
                                       const sm2_affine_t *RA,
                                       const sm2_affine_t *PA,
                                       const uint8_t *idA, size_t idA_len);
```

普通版和 Montgomery 版接口参数基本保持一致，便于在测试代码中进行差分验证和性能对比。

---

## 3. SM2 KEX 基本流程

SM2 密钥交换包含两个角色：

| 角色 | 说明 |
|---|---|
| Initiator | 发起方，通常记为 Alice |
| Responder | 响应方，通常记为 Bob |

双方各自拥有：

| 参数 | 说明 |
|---|---|
| `dA / dB` | 静态私钥 |
| `PA / PB` | 静态公钥 |
| `rA / rB` | 临时私钥 |
| `RA / RB` | 临时公钥 |
| `idA / idB` | 用户身份标识 |
| `ZA / ZB` | 与身份、公钥、曲线参数绑定的杂凑值 |
| `K` | 协商出的共享密钥 |
| `S1 / S2` | 密钥确认值 |

简化流程如下：

```text
Alice / Initiator                      Bob / Responder
---------------------------------------------------------------
生成静态密钥 dA, PA                    生成静态密钥 dB, PB
生成临时密钥 rA, RA                    生成临时密钥 rB, RB

发送 RA -----------------------------> 接收 RA
接收 RB <----------------------------- 发送 RB

计算 ZA, ZB                            计算 ZA, ZB
计算 x1', x2'                          计算 x1', x2'
计算 tA                                计算 tB
计算共享点 V                           计算共享点 U
KDF 得到共享密钥 K                     KDF 得到共享密钥 K
生成 / 验证 S1、S2                     生成 / 验证 S1、S2
```

---

## 4. 协议核心计算

SM2 密钥交换的核心目标是让 Initiator 和 Responder 在不直接传输会话密钥的情况下，基于双方的静态密钥、临时密钥和身份标识计算出一致的共享密钥。

在本项目中，SM2 KEX 的核心计算主要包括以下几个部分。

### 4.1 身份绑定值 ZA / ZB

协议首先会根据用户身份、公钥和曲线参数计算身份绑定值 `ZA` 和 `ZB`。其作用是将用户身份与密钥交换过程绑定，避免共享密钥只依赖裸公钥和临时点。

计算逻辑可以概括为：

```text
Z = SM3(ENTL || ID || a || b || Gx || Gy || Px || Py)
```

其中：

| 参数 | 含义 |
|---|---|
| `ID` | 用户身份标识 |
| `a, b` | SM2 曲线参数 |
| `Gx, Gy` | SM2 基点坐标 |
| `Px, Py` | 用户公钥坐标 |

### 4.2 x_dash 与 t 值

SM2 KEX 会从临时公钥 `RA / RB` 的 x 坐标中提取 `x_dash`，然后结合本方静态私钥和临时私钥计算标量 `t`。

项目中的计算形式可以概括为：

```text
t = d + x_dash * r mod n
```

其中：

| 符号 | 含义 |
|---|---|
| `d` | 本方静态私钥 |
| `r` | 本方临时私钥 |
| `x_dash` | 从临时公钥 x 坐标派生出的值 |
| `n` | SM2 曲线阶 |
| `t` | 后续共享点计算使用的标量 |

Initiator 侧会得到 `tA`，Responder 侧会得到 `tB`。

### 4.3 共享点与密钥派生

双方随后会计算共享椭圆曲线点。

Initiator 侧可以概括为：

```text
V = [tA](PB + [x2']RB)
```

Responder 侧可以概括为：

```text
U = [tB](PA + [x1']RA)
```

如果协议执行正确，双方最终得到的共享点应当一致。随后使用共享点坐标以及 `ZA / ZB` 进行 KDF，派生出共享密钥 `K`：

```text
K = KDF(x || y || ZA || ZB)
```

同时，协议还会计算 `S1 / S2` 等确认值，用于验证双方是否确实计算出了相同的密钥交换结果。

---

## 5. 普通路线实现

普通路线对应项目中的 Prime / 普通模 p 实现，主要文件包括：

```text
sm2_kex.c
sm2_curve.c
sm2_scalar.c
fp.c
```

该路线直接在普通有限域表示下进行计算，主要作为 baseline 和正确性参考。

普通路线的整体调用关系可以理解为：

```text
sm2_kex.c
  ↓
sm2_scalar.c
  ↓
sm2_curve.c
  ↓
fp.c
```

其中：

| 层级 | 作用 |
|---|---|
| `sm2_kex.c` | 实现普通版 SM2 密钥交换流程 |
| `sm2_scalar.c` | 实现普通域下的标量乘法 |
| `sm2_curve.c` | 实现普通域下的椭圆曲线点加、点倍等操作 |
| `fp.c` | 实现普通素域上的加减乘平方求逆等底层运算 |

在普通路线中，`RA = [rA]G`、`RB = [rB]G` 以及后续共享点计算都通过普通域标量乘法和普通域曲线点运算完成。

该路线对应的测试版本主要是：

| 版本 | 说明 |
|---|---|
| `baseline-c` | Prime 普通模 p 路线，全 C 实现 |
| `baseline-asm` | Prime 普通模 p 路线，加入 fp / fn 层 RISC-V ASM |

---

## 6. Montgomery 路线实现

Montgomery 路线是本项目的主要优化路线，主要文件包括：

```text
sm2_kex_mont.c
sm2_curve_mont.c
sm2_scalar_mont.c
fp_montgomery.c
```

该路线的协议逻辑与普通路线基本一致，但底层有限域运算改为 Montgomery 表示，以降低模乘和模约减开销。

需要注意的是，Montgomery 版 KEX 的对外接口仍然接收和输出普通 affine 点；Montgomery 表示主要用于函数内部的曲线点运算和有限域运算。这样可以保证外部接口和普通版保持一致，便于进行正确性测试和性能对比。

Montgomery 路线的整体调用关系可以理解为：

```text
sm2_kex_mont.c
  ↓
sm2_scalar_mont.c
  ↓
sm2_curve_mont.c
  ↓
fp_montgomery.c
```

其中：

| 层级 | 作用 |
|---|---|
| `sm2_kex_mont.c` | 实现 Montgomery 路线下的 SM2 密钥交换流程 |
| `sm2_scalar_mont.c` | 实现 Montgomery 域下的标量乘法 |
| `sm2_curve_mont.c` | 实现 Montgomery 域下的曲线点运算 |
| `fp_montgomery.c` | 实现 Montgomery 素域加减乘平方求逆等底层运算 |

在共享点计算中，Montgomery 路线仍然遵循：

```text
共享点 = [t](对方静态公钥 + [x_dash] 对方临时公钥)
```

区别在于中间点乘、点加和有限域运算尽量在 Montgomery 域中完成，从而减少模运算开销。

该路线对应的测试版本主要是：

| 版本 | 说明 |
|---|---|
| `mont-c` | Montgomery 路线，全 C 实现 |
| `mont-asm` | Montgomery 路线，加入 fp / fn 层 RISC-V ASM，是当前主要高性能版本 |

因此，SM2 KEX 的优化路径可以概括为：

```text
SM2 KEX
  ↓
共享点计算
  ↓
椭圆曲线标量乘法
  ↓
曲线点加 / 点倍
  ↓
有限域 Montgomery 运算
  ↓
RISC-V ASM 优化
```

---

## 7. 正确性验证

当前正确性测试使用 `src/test_correct.c`，测试命令为：

```bash
cd src
make clean
make correctness
```

`make correctness` 会分别运行：

```text
test_correct_prime_c
test_correct_prime_fp_asm
test_correct_mont_c
test_correct_mont_asm
```

对应四个版本：

| 版本 | 说明 |
|---|---|
| `baseline-c` | Prime 普通模 p 路线，全 C 实现 |
| `baseline-asm` | Prime 普通模 p 路线，fp / fn ASM 实现 |
| `mont-c` | Montgomery 路线，全 C 实现 |
| `mont-asm` | Montgomery 路线，fp / fn ASM 实现 |

测试内容包括固定 SM2 KEX 测试向量、Alice / Bob 共享密钥一致性、`S1 / S2` 校验值验证，以及 10 轮随机 SM2 KEX 测试。

当前四个版本均已通过测试，说明 Montgomery 优化和 RISC-V ASM 优化没有破坏 SM2 密钥交换协议正确性。

---

## 8. KEX 性能结果摘要

当前性能测试主要比较四个版本：`baseline-c`、`baseline-asm`、`mont-c`、`mont-asm`。

其中，`mont-asm` 是当前最快实现路径。相比初始 `baseline-c`，整体运行时间从 12.13s 降低到 5.97s，约获得 **2.03x 加速**。

在 KEX 关键路径上，`KEX Initiator verify` 从 6,333,281 cycles 降低到 3,165,830 cycles，接近 **2.00x 加速**。

这说明 Montgomery 路线是主要优化来源，进一步结合 RISC-V ASM 后，优化效果能够传导到 SM2 KEX 协议层。

完整性能测试结果见 `docs/benchmark.md`。

---

## 9. KEX 专项 perf workload

除 `make perf-all` 外，项目还提供 `test_perf.c`，用于给 Linux `perf stat` / `perf record` 提供更干净的 KEX workload。它支持 `genra`、`scalar`、`responder`、`initiator`、`full` 和 `all` 等模式。

常用命令包括：

```bash
make perf-test-perf-genra
make perf-test-perf-scalar
make perf-test-perf-responder
make perf-test-perf-initiator
```

更详细的 perf 使用说明见 `docs/benchmark.md`。

---

## 10. 与 SM4 安全通信 demo 的关系

SM2 KEX 的输出是双方协商得到的共享密钥 `K`。在本项目中，SM4 安全通信 demo 已经基于该共享密钥实现了上层消息加密演示。

整体流程如下：

```text
SM2 KEX
  ↓
双方协商得到 shared key
  ↓
KDF 派生 SM4 会话密钥
  ↓
使用 SM4-CBC 加密用户输入消息
  ↓
使用相同会话密钥解密密文
  ↓
验证解密结果是否等于原始明文
```

该 demo 的作用是展示优化后的 SM2 密钥交换结果如何用于实际安全通信场景。需要注意的是，SM4 demo 不直接调用底层 RISC-V 汇编函数，也不重新实现 SM2，而是调用上层 SM2 KEX 接口获取 shared key，再完成密钥派生和 SM4-CBC 加解密。

CLI 菜单中对应功能为：

```text
2. Run SM4 secure message demo
```

因此，SM2 KEX 是项目底层优化和上层安全通信演示之间的连接点。

---

## 11. 小结

SM2 KEX 是本项目连接底层优化和协议层应用的关键模块。

本项目通过：

- Prime 普通模 p 路线与 Montgomery 路线对比；
- C 实现与 RISC-V ASM 实现对比；
- `baseline-c`、`baseline-asm`、`mont-c`、`mont-asm` 四版本正确性测试；
- SM2 KEX 固定测试向量与随机测试；
- Linux `perf` 性能评测；
- SM4-CBC 安全通信 demo；

验证了底层 Montgomery 运算、标量乘法和 RISC-V ASM 优化能够传导到 SM2 密钥交换协议层，并进一步支撑上层 SM4 安全通信演示。

当前结果表明，`mont-asm` 是项目中最快的 SM2 KEX 实现路径，相比初始 `baseline-c` 在整体 perf 测试中约获得 **2.03x 加速**。
