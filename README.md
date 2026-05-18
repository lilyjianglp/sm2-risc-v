# RV-GM-Secure：SM2 优化后端与 SM4 安全通信演示

本分支在原有 RISC-V SM2 优化项目基础上，新增了一个完整的安全通信演示层。

项目核心仍然是 **RISC-V 平台上的 SM2 密钥交换优化**。新增的 SM4 安全通信模块不会重新实现 SM2，也不会修改底层汇编优化逻辑，而是调用上层 SM2 密钥交换结果，将双方协商得到的 shared key 派生为 SM4 会话密钥，并完成消息加密、解密和一致性验证。

## 本分支新增内容

本分支主要新增两个部分：

```text
1. CLI 安全通信演示
2. React 前端展示页面
```

CLI 是真实后端运行入口，用于在 RISC-V 板端执行 SM2 密钥交换、SM4 加解密、正确性测试和性能评测。

前端是展示层，用于可视化说明 SM2 优化链路、SM4 安全信道、板端实测数据和运行命令。

## 核心密码链路

完整链路如下：

```text
SM2 密钥交换
-> Montgomery 有限域运算后端
-> Pornin 全模逆 RISC-V 汇编优化
-> 生成双方一致的 shared key，验证 KA == KB
-> SM3/KDF 派生 128-bit SM4 会话密钥
-> SM4-CBC 加密和解密用户消息
```

其中：

- **SM2** 是本项目的核心优化对象；
- **Montgomery 后端** 用于加速有限域模乘、模平方和模约减；
- **Pornin 全模逆汇编** 用于优化 `fp_mont_inv` 关键路径；
- **SM4** 是基于 SM2 shared key 的上层安全通信演示；
- **SM3/KDF** 用于从 SM2 shared key 派生 SM4 会话密钥。

## CLI 后端演示

CLI 代码位置：

```text
sm2-risc-v-main/src/
```

主要文件：

```text
sm2-risc-v-main/src/sm4_channel.h
sm2-risc-v-main/src/sm4_channel.cpp
sm2-risc-v-main/src/demo/secure_message_demo.cpp
sm2-risc-v-main/src/demo/sm4_demo.md
sm2-risc-v-main/src/Makefile
```

编译和运行：

```bash
cd sm2-risc-v-main/src
make secure_message_demo
./secure_message_demo
```

CLI 菜单：

```text
=========================================
 RV-GM-Secure: RISC-V GM Crypto Backend
=========================================
1. Run SM2 key exchange demo
2. Run SM4 secure message demo
3. Run correctness test
4. Run performance benchmark
0. Exit
=========================================
Please select:
```

### 1. SM2 密钥交换演示

运行优化后的 SM2 密钥交换流程，调用 Montgomery 后端和 Pornin 全模逆汇编优化路径，生成双方 shared key，并验证双方密钥一致：

```text
[SM2 KEX] Running optimized SM2 key exchange...
[SM2 KEX] Shared key generated successfully.
[SM2 KEX] KA == KB: true
```

### 2. SM4 安全通信演示

用户输入一段明文消息后，程序执行以下流程：

```text
输入明文
-> 运行 optimized SM2 key exchange
-> 检查 KA == KB
-> 使用 SM3/KDF 派生 SM4 会话密钥
-> 生成 16 字节 IV
-> 使用 SM4-CBC 加密消息
-> 使用同一 key 和 IV 解密密文
-> 验证解密结果是否等于原始明文
-> 输出 PASS 或 FAIL
```

运行示例：

```text
Please input message:
hello im Aria

[SM2 KEX] Running optimized SM2 key exchange...
[SM2 KEX] Shared key generated successfully.
[SM2 KEX] KA == KB: true

[KDF] Deriving SM4 session key from SM2 shared key...
[SM4] IV: ...
[SM4] Ciphertext: ...

[SM4] Decrypting message...
[SM4] Plaintext: hello im Aria

[PASS] Secure message demo success.
```

### 3. 正确性测试

运行 Pornin 全模逆汇编相关正确性测试，用于验证汇编实现和 C 参考实现的一致性：

```bash
make run-pornin-inv-asm USE_FULL_INV_ASM=1
```

板端验证结果示例：

```text
[OK] inner31 asm/c reference tests passed
[OK] fixed tests passed
[OK] random tests passed
[OK] zero test passed
all tests passed
```

### 4. 性能评测

运行板端性能评测：

```bash
make perf-pornin-asm USE_FULL_INV_ASM=1
```

RISC-V Muse Pi Pro 板端实测结果示例：

```text
fp_mont_mul (ASM)                      : mean=330 cycles
fp_mont_sqr (ASM)                      : mean=298 cycles
fp_mont_inv (Pornin full ASM)          : mean=17907 cycles
KEX Responder Mont-C (full)            : mean=3128943 cycles
KEX Initiator Mont-C (verify)          : mean=2846039 cycles
IPC                                    : 1.36 insn/cycle
```

## 前端展示页面

前端代码位置：

```text
frontend/
```

主要文件：

```text
frontend/index.html
frontend/package.json
frontend/package-lock.json
frontend/vite.config.js
frontend/README.md
frontend/src/main.jsx
frontend/src/styles.css
frontend/src/demoData.js
```

本地运行：

```bash
cd frontend
npm install
npm run dev
```

Windows PowerShell 下建议使用：

```powershell
cd frontend
npm.cmd install
npm.cmd run dev -- --host 127.0.0.1 --port 3000
```

浏览器打开：

```text
http://127.0.0.1:3000/
```

前端包含 6 个页面：

```text
1. SM2 优化概览
2. SM2 密钥交换
3. SM4 安全信道
4. RISC-V 后端链路
5. 性能基准评测
6. 运行指南
```

## 前端与 CLI 的关系

前端不是密码算法实现，也不会直接调用 RISC-V 板端程序。

前端的作用是展示和解释：

- SM2 密钥交换优化链路；
- Montgomery 有限域后端；
- Pornin 全模逆 RISC-V 汇编；
- SM3/KDF 派生 SM4 会话密钥；
- SM4-CBC 安全通信流程；
- RISC-V 板端性能测试结果；
- CLI 编译和运行命令。

前端展示数据来自板端 CLI 实际运行结果，保存在：

```text
frontend/src/demoData.js
```

如果用户输入内容匹配已有板端样例，前端展示真实板端数据。

如果用户输入新内容，前端进入预览模式，仅用于展示交互效果。

真实密码后端仍以 CLI 为准：

```text
sm2-risc-v-main/src/secure_message_demo
```

## 分层设计说明

本分支保持了清晰的软件分层：

```text
前端 / CLI 展示层
-> SM2 密钥交换协议层
-> SM2 标量乘法层
-> Montgomery 有限域后端
-> Pornin 全模逆 RISC-V 汇编
-> SM3/KDF 密钥派生层
-> SM4-CBC 安全通信层
```

SM4 安全通信模块只调用上层 SM2 密钥交换结果，不修改底层 RISC-V 汇编优化代码，也不重新实现 SM2。

## 新增或修改文件

```text
README.md
frontend/
sm2-risc-v-main/src/sm4_channel.h
sm2-risc-v-main/src/sm4_channel.cpp
sm2-risc-v-main/src/demo/
sm2-risc-v-main/src/Makefile
```

## 总结

本分支将原有 SM2 RISC-V 优化项目扩展为一个可演示的国密安全通信系统：

```text
optimized SM2 key exchange
+ Montgomery / Pornin RISC-V optimized backend
+ SM3/KDF key derivation
+ SM4-CBC secure message demo
+ CLI interaction
+ frontend visualization
```

核心贡献可以概括为：

```text
以 RISC-V 平台 SM2 密钥交换优化为核心，
在上层构建 SM4 安全通信演示，
并提供 CLI 与前端两种展示入口。
```
