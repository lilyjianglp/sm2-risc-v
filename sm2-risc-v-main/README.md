# SM2-RISC-V: High-Performance Assembly Optimization

This project presents a high-performance implementation of the SM2 (Chinese National Standard) elliptic curve cryptography algorithm optimized for the RISC-V 64-bit (RV64) architecture.

## Overview

Elliptic Curve Cryptography (ECC) performance is largely dominated by prime field arithmetic operations such as modular addition, subtraction, multiplication, and reduction.

This project accelerates these critical operations by replacing standard C implementations with hand-optimized RV64 assembly, significantly improving execution efficiency on RISC-V processors.

These optimized arithmetic routines are integrated into the SM2 key exchange protocol, improving the efficiency of elliptic curve computations during key negotiation.


## 🛠️ Getting Started

### Prerequisites

To build and run the project, you need the following:

- **RISC-V Cross-Compilation Toolchain**
- **QEMU** for RISC-V simulation

### Installation

On Ubuntu/Debian, install the necessary dependencies as follows:

```bash
sudo apt install gcc-riscv64-linux-gnu qemu-user
```

### Compilation

The project offers two versions for building:

**1. Assembly Optimized Version** — uses hand-written assembly code for optimized performance:

```bash
make USE_FP_ASM=1
```

**2. Pure C Reference Version** — uses the standard C implementation:

```bash
make USE_FP_ASM=0
```

### Execution

To run the benchmark suite and test the optimizations:

```bash
make run
```

---

## 🧠 Technical Implementation Details

### 1. Branchless Carry Chain Optimization

In cryptographic operations, branches like `if (carry)` can cause mispredictions. Instead, we use the `sltu` (Set Less Than Unsigned) instruction to capture the carry bit in a branchless manner:

```asm
add   a0, a0, t1      # Add partial product
sltu  t2, a0, t1      # Capture carry bit into t2 in 1 cycle
```

### 2. Squaring Optimization

The squaring operation, a critical part of ECC, leverages symmetry: $a_i \cdot a_j = a_j \cdot a_i$. This symmetry allows us to compute mixed products once and double them, reducing the number of multiplication instructions from **16 to 10** for a 256-bit operation.

```asm
# Optimized squaring for 256-bit prime field
```

---
## Benchmark Results
Benchmarks were measured using the `rdcycle` instruction under QEMU RISC-V emulator.
Each test was executed multiple times and averaged to obtain stable results.

All benchmarks were executed using **QEMU RISC-V simulation** with cycle counting.  
The results compare three implementations:

- **C** — reference C implementation  
- **RV ASM** — hand-optimized RISC-V assembly implementation  
- **Mont** — Montgomery optimized version

---

### Basic Arithmetic Operations (cycles)

| Operation   | C       | RV ASM | Montgomery |
|-------------|---------|--------|-----------|
| fp_reduce   | 241.25  | 217.65 | 213.05 |
| fp_mul      | 341.90  | 208.75 | 204.85 |
| fp_sqr      | 261.75  | 192.00 | 185.00 |

---

### High-Level ECC Operations (cycles)

| Operation   | C        | RV ASM     | Montgomery |
|-------------|----------|-----------|-----------|
| ScalarMul   | 2,099,501.40 | 1,618,641.90 | 1,526,444.10 |
| KEX Verify  | 4,486,842.20 | 3,473,295.25 | 3,332,707.30 |

---

### Speedup (vs C)

| Operation | RV ASM | Montgomery |
|----------|--------|-----------|
| fp_mul | **1.63×** | **1.66×** |
| ScalarMul | **1.29×** | **1.37×** |
| KEX Verify | **1.29×** | **1.34×** |
## 🔧 Project Structure

```
/sm2-riscv/
├── README.md                      # Project documentation
├── src/                           # Source code directory
│   ├── asm/                       # Assembly optimizations
│   │   ├── fp_add.S               # Constant-time addition
│   │   ├── fp_sub.S               # Constant-time subtraction
│   │   ├── fp_mul.S               # Optimized multiplication
│   │   ├── fp_sqr.S               # Optimized squaring
│   │   ├── fp_neg.S               # Constant-time negation
│   │   ├── fp_reduce.S            # Field reduction
│   │   ├── fn_mul_montgomery.s    # Montgomery multiplication (ASM)
│   │   └── fn_reduce_Montgomery.s # Montgomery reduction (ASM)
│   ├── bench.sh                   # Benchmark script
│   ├── fn.c                       # Scalar field function implementations
│   ├── fn.h                       # Header for scalar field functions
│   ├── fn_Montgomery.c            # Montgomery multiplication in C
│   ├── fp.c                       # Prime field operations
│   ├── fp.h                       # Header for prime field
│   ├── Makefile                   # Build system
│   ├── read_cycle.s               # Cycle reading for performance analysis
│   ├── sm2_curve.c                # Elliptic curve point operations
│   ├── sm2_curve.h                # Header for elliptic curve
│   ├── sm2_kex.c                  # Key exchange protocol
│   ├── sm2_kex.h                  # Header for key exchange
│   ├── sm2_scalar.c               # Scalar operations
│   ├── sm2_scalar.h               # Header for scalar operations
│   └── test.c                     # Benchmarking and validation tests

```

---

## Authors

- Liangping Jiang
- Yiwen Wang

### Research Direction

- RISC-V Architecture
- Cryptographic Engineering
- SM2 Implementation
