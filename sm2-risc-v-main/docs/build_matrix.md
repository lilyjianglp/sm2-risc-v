# 构建版本说明

本项目为了对比不同实现路径下 SM2 密钥交换的正确性和性能，保留了四个构建版本：

- `baseline-c`
- `baseline-asm`
- `mont-c`
- `mont-asm`

这四个版本用于分别观察两类优化效果：

1. Prime 普通模 p 路线到 Montgomery 路线的算法优化效果；
2. C 实现到 RISC-V ASM 实现的底层优化效果。

---

## 1. 构建版本总览

| 版本 | 算法路线 | 实现方式 | 是否使用汇编 | 说明 |
|---|---|---|---|---|
| `baseline-c` | Prime 普通模 p 路线 | fp C + fn C | 否 | 初版全 C baseline |
| `baseline-asm` | Prime 普通模 p 路线 | fp ASM + fn ASM | 是 | 初版路线加入 RISC-V 汇编优化 |
| `mont-c` | Montgomery 路线 | fp C + fn C | 否 | Montgomery 优化路线的 C 版本 |
| `mont-asm` | Montgomery 路线 | fp ASM + fn ASM | 是 | Montgomery + RISC-V ASM 最终优化版本 |

其中：

- `baseline-c` 和 `mont-c` 不编译任何 `.S` 汇编文件；
- `baseline-asm` 和 `mont-asm` 会编译对应的 RISC-V 汇编优化文件；
- `mont-asm` 是当前项目的最终高性能版本。

---

## 2. baseline-c

`baseline-c` 是初版 Prime 普通模 p 路线的全 C 实现。

### 编译内容

```text
fp.c
sm2_curve.c
sm2_scalar.c
sm2_kex.c
fn_Montgomery.c
```

### 是否使用汇编

不使用任何 `.S` 汇编文件。

### 运行命令

```bash
make baseline-c
```

### 作用

`baseline-c` 作为项目的初始基准版本，用于：

- 正确性参考；
- 性能对比基准；
- 计算后续优化版本的加速比。

---

## 3. baseline-asm

`baseline-asm` 是 Prime 普通模 p 路线加入 RISC-V ASM 优化后的版本。

### 编译内容

```text
fp.c
sm2_curve.c
sm2_scalar.c
sm2_kex.c
fn_Montgomery.c
asm/fn_mont.S
asm/fp_add.S
asm/fp_sub.S
asm/fp_neg.S
asm/fp_mul.S
asm/fp_sqr.S
asm/fp_reduce.S
```

### 是否使用汇编

使用 RISC-V 汇编文件。

### 运行命令

```bash
make baseline-asm
```

### 作用

`baseline-asm` 用于观察在 Prime 普通模 p 路线下加入 RISC-V ASM 后，对整体 SM2 密钥交换性能是否有提升。

---

## 4. mont-c

`mont-c` 是 Montgomery 路线的全 C 实现。

### 编译内容

```text
fp_montgomery.c
sm2_curve_mont.c
sm2_scalar_mont.c
sm2_kex_mont.c
fn_Montgomery.c
```

### 是否使用汇编

不使用任何 `.S` 汇编文件。

### 运行命令

```bash
make mont-c
```

### 作用

`mont-c` 用于单独观察 Montgomery 算法路线本身带来的性能提升。

---

## 5. mont-asm

`mont-asm` 是 Montgomery 路线加入 RISC-V ASM 优化后的版本，也是当前项目的最终高性能版本。

### 编译内容

```text
fp_montgomery.c
sm2_curve_mont.c
sm2_scalar_mont.c
sm2_kex_mont.c
fn_Montgomery.c
asm/fn_mont.S
asm/fp_mont_add.S
asm/fp_mont_sub.S
asm/fp_mont_neg.S
asm/fp_mont_mul.S
asm/fp_mont_sqr.S
asm/fp_mont_to_from.S
```

### 是否使用汇编

使用 RISC-V 汇编文件。

### 运行命令

```bash
make mont-asm
```

### 作用

`mont-asm` 是最终优化版本，用于展示 Montgomery 表示与 RISC-V 汇编优化结合后的性能效果。

---

## 6. 常用命令

### 编译并运行单个版本

```bash
make baseline-c
make baseline-asm
make mont-c
make mont-asm
```

### 编译全部版本

```bash
make clean
make all
```

### 运行正确性测试

```bash
make correctness
```

### 运行性能测试

```bash
make perf-all
```

---

## 7. 对比关系

四个版本形成如下对比关系：

| 对比关系 | 说明 |
|---|---|
| `baseline-c` vs `baseline-asm` | 对比 Prime 路线下 ASM 优化效果 |
| `baseline-c` vs `mont-c` | 对比 Montgomery 算法路线优化效果 |
| `mont-c` vs `mont-asm` | 对比 Montgomery 路线下 ASM 优化效果 |
| `baseline-c` vs `mont-asm` | 对比初版 baseline 与最终优化版本的整体差异 |

---

## 8. 当前结论

当前测试结果表明：

- 四个版本均可成功构建和运行；
- `baseline-c` 和 `mont-c` 是不含汇编的 C 版本；
- `baseline-asm` 和 `mont-asm` 是带 RISC-V 汇编优化的版本；
- `mont-asm` 是当前最快版本，也是最终展示的优化版本。

