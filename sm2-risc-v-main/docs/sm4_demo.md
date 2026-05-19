# SM4 安全通信演示模块

本文档说明面向 64 位 RISC-V 平台的 SM2 密钥交换算法高效实现项目中的 SM4 安全通信演示模块，包括模块作用、代码结构、运行方式、交互流程和分层设计。

SM4 demo 的作用是展示：优化后的 SM2 密钥交换结果如何进一步用于上层安全通信。程序首先运行 SM2 密钥交换，得到通信双方一致的共享密钥，然后通过 KDF 派生 SM4 会话密钥，最后使用 SM4-CBC 对用户输入的消息进行加密和解密验证。

---

## 1. 模块定位

SM4 安全通信演示模块位于项目上层，不直接修改底层 RISC-V 汇编优化代码。

整体调用关系如下：

```text
底层优化层：
Montgomery 有限域运算 / 标量乘法 / Pornin 模逆 / RISC-V ASM
        ↓
SM2 密钥交换层：
生成双方一致 shared key
        ↓
KDF 派生层：
从 SM2 shared key 派生 SM4 会话密钥
        ↓
SM4 安全通信层：
SM4-CBC 加密 / 解密用户消息
        ↓
CLI 交互层：
用户输入消息，展示加密、解密和验证结果
```

因此，SM4 demo 的核心不是重新优化 SM4，也不是重新实现 SM2，而是把已经完成的 SM2 KEX 结果用于一个实际的安全通信场景。

---

## 2. 代码结构

SM4 demo 相关文件主要包括：

```text
src/
├── sm4_channel.cpp
├── sm4_channel.h
└── demo/
    └── secure_message_demo.cpp
```

其中：

| 文件 | 说明 |
|---|---|
| `sm4_channel.h` | 定义 SM2 KEX 结果结构、SM4 安全通信结果结构，以及对外接口 |
| `sm4_channel.cpp` | 实现 SM2 KEX 封装、KDF、IV 生成、SM4-CBC 加密和解密 |
| `demo/secure_message_demo.cpp` | CLI 菜单入口，负责用户交互和结果打印 |

---

## 3. 对外接口

`sm4_channel.h` 中提供了以下主要结构和接口。

### 3.1 结果结构体

```cpp
struct Sm2KexResult {
    bool success = false;
    std::vector<uint8_t> shared_key_a;
    std::vector<uint8_t> shared_key_b;
};

struct SecureMessageResult {
    bool kex_success = false;
    bool decrypt_success = false;

    std::string plaintext;
    std::string decrypted_text;

    std::vector<uint8_t> shared_key_digest;
    std::vector<uint8_t> sm4_key_digest;
    std::vector<uint8_t> iv;
    std::vector<uint8_t> ciphertext;
};
```

### 3.2 核心接口

```cpp
Sm2KexResult run_sm2_key_exchange_demo(bool use_optimized_sm2);

std::vector<uint8_t> derive_sm4_key(
    const std::vector<uint8_t>& shared_key
);

std::vector<uint8_t> generate_random_iv();

std::vector<uint8_t> sm4_encrypt_message(
    const std::vector<uint8_t>& key,
    const std::vector<uint8_t>& iv,
    const std::string& plaintext
);

std::string sm4_decrypt_message(
    const std::vector<uint8_t>& key,
    const std::vector<uint8_t>& iv,
    const std::vector<uint8_t>& ciphertext
);

SecureMessageResult run_sm4_secure_message_demo(
    const std::string& message,
    bool use_optimized_sm2
);

std::string bytes_to_hex(const std::vector<uint8_t>& bytes);
```

需要注意的是，当前实现中 `run_sm2_key_exchange_demo(bool use_optimized_sm2)` 的参数没有在运行时切换不同后端；实际使用普通路径还是 Montgomery 优化路径，由编译时宏和 Makefile 链接内容决定。

---

## 4. 功能流程

SM4 secure message demo 的主流程如下：

```text
用户输入明文消息
        ↓
运行 SM2 key exchange
        ↓
生成 shared key
        ↓
验证 KA == KB
        ↓
通过 KDF 派生 SM4 session key
        ↓
生成 IV
        ↓
使用 SM4-CBC 加密消息
        ↓
使用相同 key 和 IV 解密密文
        ↓
比较解密结果与原始明文
        ↓
输出 PASS / FAIL
```

其中：

| 步骤 | 说明 |
|---|---|
| SM2 KEX | 负责让通信双方生成一致的 shared key |
| KDF | 从 shared key 派生 16 字节 SM4 会话密钥 |
| IV | 为 SM4-CBC 生成 16 字节初始化向量 |
| SM4-CBC | 对用户输入消息进行加密和解密 |
| PASS / FAIL | 判断解密结果是否与原始输入一致 |

---

## 5. SM2 KEX 封装

`sm4_channel.cpp` 中通过 `run_sm2_key_exchange_demo()` 封装 SM2 密钥交换流程。

该函数内部会：

1. 使用固定测试私钥和临时私钥；
2. 生成 Alice / Bob 的静态公钥；
3. 生成双方临时公钥 `RA / RB`；
4. 分别运行 Initiator 和 Responder 侧 KEX；
5. 检查双方 shared key 是否一致；
6. 检查 `S1 / S2` 确认值是否一致；
7. 返回 `shared_key_a` 和 `shared_key_b`。

SM4 demo 在进入加密阶段前，会检查：

```text
kex.success == true
shared_key_a.size() == 32
shared_key_a == shared_key_b
```

只有上述条件满足时，才继续派生 SM4 会话密钥。

---

## 6. SM4 会话密钥派生

SM4 的密钥长度为 16 字节。当前实现不会直接把 SM2 shared key 截断作为 SM4 key，而是使用带标签的 KDF 方式派生。

派生方式为：

```text
SM4_key = SM3(shared_key || "RV-GM-Secure-SM4-Demo")[0:16]
```

其中：

| 字段 | 说明 |
|---|---|
| `shared_key` | SM2 KEX 输出的共享密钥 |
| `"RV-GM-Secure-SM4-Demo"` | demo 专用 KDF 标签 |
| `SM3(...)` | 对拼接后的输入做 SM3 哈希 |
| `[0:16]` | 取 SM3 输出前 16 字节作为 SM4 key |

这样可以避免直接使用 shared key 作为对称加密密钥，也方便区分 demo 场景下的密钥用途。

---

## 7. IV 生成

SM4-CBC 需要 16 字节 IV。当前实现会优先从：

```text
/dev/urandom
```

读取随机字节作为 IV。若读取失败，则使用基于时间和地址扰动的 fallback 方式生成 IV。

该 fallback 主要用于 demo 容错；正式安全场景中应优先依赖系统安全随机数。

---

## 8. SM4-CBC 加密与解密

当前实现使用 SM4-CBC 模式处理用户输入消息。

加密流程为：

```text
明文字符串
  ↓
PKCS#7 padding
  ↓
每个明文块与上一轮密文块 / IV 异或
  ↓
SM4 block encryption
  ↓
输出 ciphertext
```

解密流程为：

```text
ciphertext
  ↓
SM4 block decryption
  ↓
与上一轮密文块 / IV 异或
  ↓
去除 PKCS#7 padding
  ↓
恢复明文字符串
```

`run_sm4_secure_message_demo()` 会在完成解密后比较：

```text
decrypted_text == message
```

如果一致，则设置：

```text
decrypt_success = true
```

---

## 9. 编译运行

在项目的 `src` 目录下执行：

```bash
cd src
make secure_message_demo
./secure_message_demo
```

也可以直接使用 Makefile 中的运行目标：

```bash
make run-secure-message-demo
```

当前 `secure_message_demo` 目标链接的是 Montgomery + fp ASM + Pornin full inversion ASM 的优化路径，包括：

```text
fp_montgomery.o
fp_pornin_inv.o
sm2_curve_mont.o
sm2_scalar_mont.o
sm2_kex_mont.o
fp_mont_*.o
pornin_full_inv.o
fn_mont.o
```

---

## 10. CLI 菜单说明

运行 `./secure_message_demo` 后，会进入如下菜单：

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
```

各选项含义如下：

| 选项 | 功能 |
|---|---|
| `1` | 运行 SM2 密钥交换 demo，验证双方 shared key 是否一致 |
| `2` | 运行 SM4 安全通信 demo，对用户输入消息进行加密和解密 |
| `3` | 运行当前优化路径下的 Pornin full inversion ASM 正确性测试 |
| `4` | 运行当前优化路径下的 Pornin ASM 性能测试 |
| `0` | 退出程序 |

SM4 安全通信演示对应菜单选项：

```text
2. Run SM4 secure message demo
```

---

## 11. 运行示例

示例输入：

```text
hello im Aria
```

运行输出示例：

```text
Please input message:
hello im Aria

[SM2 KEX] Running optimized SM2 key exchange...
[SM2 KEX] Shared key generated successfully.
[SM2 KEX] KA == KB: true

[KDF] Deriving SM4 session key from SM2 shared key...
[KDF] SM4 key: 3b5342a6ab9c60d09486788aae0192ba
[SM4] IV: e2dc8e2d10712bdcaa6ed9fefc5558
[SM4] Encrypting message...
[SM4] Ciphertext: ef7a092e929e4f0256abb22a2ac5a1cc

[SM4] Decrypting message...
[SM4] Plaintext: hello im Aria

[PASS] Secure message demo success.
```

该结果说明：

1. SM2 KEX 成功运行；
2. Alice 和 Bob 得到的共享密钥一致，即 `KA == KB: true`；
3. 程序成功从 shared key 派生 SM4 会话密钥；
4. SM4-CBC 成功完成消息加密；
5. 解密后的明文与用户输入一致；
6. 安全通信演示通过。

---

## 12. 分层设计说明

本模块保持了清晰的分层结构：

| 层级 | 作用 |
|---|---|
| RISC-V ASM / Montgomery 层 | 提供底层高性能有限域运算 |
| SM2 KEX 层 | 基于底层运算生成双方一致 shared key |
| KDF 层 | 从 shared key 派生 SM4 session key |
| SM4-CBC 层 | 加密和解密用户消息 |
| CLI 层 | 提供用户输入和结果展示 |

这种设计的好处是：

1. SM4 demo 不破坏底层优化代码；
2. SM2 KEX 可以独立测试和 benchmark；
3. SM4 安全通信流程可以作为上层应用演示；
4. CLI 输出直观，适合项目展示和答辩。

---

## 13. 注意事项

当前 demo 会在终端中输出 SM4 key、IV 和 ciphertext，主要是为了便于展示和调试。在真实应用场景中，会话密钥不应直接打印到终端。

此外，SM4 demo 的重点是展示完整安全通信链路：

```text
SM2 密钥交换 → KDF 派生密钥 → SM4-CBC 加密通信
```

它不是 SM4 性能优化模块，也不承担底层 RISC-V 汇编优化任务。

---

## 14. 小结

SM4 安全通信 demo 是 SM2-RISC-V 项目的上层应用演示模块。

它通过：

1. 运行优化后的 SM2 密钥交换；
2. 验证 `KA == KB`；
3. 从 shared key 派生 SM4 session key；
4. 使用 SM4-CBC 加密和解密用户消息；
5. 验证解密结果与原始明文一致；

展示了本项目从底层 SM2 优化到上层安全通信应用的完整链路。