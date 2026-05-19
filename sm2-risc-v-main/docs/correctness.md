# 正确性测试说明

本项目新增了独立的 SM2 正确性测试流程，用于验证四个构建版本在 SM2 密钥交换协议层面的输出是否正确。

当前正确性测试覆盖四个版本：

- `baseline-c`
- `baseline-asm`
- `mont-c`
- `mont-asm`

---

## 1. 测试目标

正确性测试的目标是验证：

1. 四个版本都能成功完成 SM2 密钥交换；
2. Alice 和 Bob 双方能够得到一致的共享密钥；
3. S1 / S2 确认值能够正确生成和验证；
4. 固定测试向量下，四个版本输出完全一致；
5. 随机测试下，四个版本均能稳定通过。

---

## 2. 测试命令

在 `src/` 目录下执行：

```bash
make clean
make correctness
```

也可以在项目根目录下使用一键脚本：

```bash
bash scripts/run_correctness.sh
```

---

## 3. 测试程序

`make correctness` 会依次构建并运行四个正确性测试程序：

| 测试程序 | 对应版本 | 说明 |
|---|---|---|
| `test_correct_prime_c` | `baseline-c` | Prime 路线，fp C + fn C |
| `test_correct_prime_fp_asm` | `baseline-asm` | Prime 路线，fp ASM + fn ASM |
| `test_correct_mont_c` | `mont-c` | Montgomery 路线，fp C + fn C |
| `test_correct_mont_asm` | `mont-asm` | Montgomery 路线，fp ASM + fn ASM |

每个测试程序只链接自身对应路线的源码，避免把 Prime 路线和 Montgomery 路线同时链接到一个可执行文件中导致符号冲突。

---

## 4. 测试内容

正确性测试主要包括两部分：

### 4.1 固定 SM2 KEX 测试

固定测试使用固定的：

- Alice 私钥 `dA`
- Bob 私钥 `dB`
- Alice 临时随机数 `rA`
- Bob 临时随机数 `rB`
- Alice ID
- Bob ID

测试流程包括：

1. 生成 Alice 和 Bob 的长期公钥；
2. 生成 Alice 和 Bob 的临时公钥；
3. Alice 执行第一阶段密钥计算；
4. Bob 执行响应方密钥计算；
5. Alice 执行第二阶段确认；
6. 检查双方共享密钥是否一致；
7. 检查 S1 / S2 确认值是否一致。

### 4.2 随机 SM2 KEX 测试

每个版本还会执行 10 轮随机 SM2 KEX 测试，用于验证不同输入下协议流程是否稳定。

每一轮随机测试都会重新生成：

- Alice 私钥；
- Bob 私钥；
- Alice 临时随机数；
- Bob 临时随机数。

---

## 5. 当前测试结果

四个版本均通过固定 SM2 KEX 测试和 10 轮随机 SM2 KEX 测试。

| 版本 | 固定 KEX 测试 | 随机 KEX 测试 | 结果 |
|---|---|---|---|
| `baseline-c` | PASS | 10/10 PASS | PASS |
| `baseline-asm` | PASS | 10/10 PASS | PASS |
| `mont-c` | PASS | 10/10 PASS | PASS |
| `mont-asm` | PASS | 10/10 PASS | PASS |

固定测试向量下，四个版本输出完全一致：

```text
SharedKey = db0e9fdca7c9b8cde2edf89e6ec184ffeb3c34dbe50bc8ad72c11f7735f4108f
S1        = 53752ebb074279b5c909d88c4ae9bb611e5dc90dba3eb079ce1133d62ff489d3
S2        = 431f7b64ceaa65f8cadffcc5b6056530c2139b4f2fd26b2fb6d05e251c709b85
```

---

## 6. 输出示例

运行 `make correctness` 后，可以看到类似输出：

```text
=== correctness: baseline-c / Prime / fp C / fn C ===
[PASS] SM2 correctness test passed

=== correctness: baseline-asm / Prime / fp ASM / fn ASM ===
[PASS] SM2 correctness test passed

=== correctness: mont-c / Montgomery / fp C / fn C ===
[PASS] SM2 correctness test passed

=== correctness: mont-asm / Montgomery / fp ASM / fn ASM ===
[PASS] SM2 correctness test passed

[PASS] correctness workflow completed
```

---

## 7. 设计说明

一开始尝试过在同一个测试程序中同时链接 Prime 路线和 Montgomery 路线进行直接对拍，但由于两套路线上存在部分同名函数，会导致链接阶段出现 `multiple definition` 错误。

因此最终采用当前设计：

- 四个版本分别编译独立的正确性测试程序；
- 每个程序只链接自身实现路径；
- 通过固定测试向量输出的 `SharedKey`、`S1`、`S2` 来检查四个版本协议输出是否一致。

这种方式可以避免符号冲突，同时仍然能够验证四个实现版本在 SM2 KEX 协议层面的正确性。

---

## 8. 结论

正确性测试结果表明：

- `baseline-c`、`baseline-asm`、`mont-c`、`mont-asm` 四个版本均能正确完成 SM2 密钥交换；
- 四个版本在固定测试向量下得到完全一致的 `SharedKey`、`S1` 和 `S2`；
- Montgomery 优化没有破坏协议正确性；
- RISC-V ASM 优化没有破坏协议正确性；
- 当前项目已经具备可复现的 SM2 正确性测试流程。
