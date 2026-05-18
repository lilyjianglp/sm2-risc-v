# SM4 安全通信演示模块

这个文件夹集中放置 SM4 安全通信演示相关的“展示入口”和说明文档，方便在板子的 `src/` 目录里提交和查找。

## 文件说明

- `secure_message_demo.cpp`: CLI 交互演示入口，包含菜单和用户输入流程。
- `sm4_demo.md`: 模块设计、SM2/SM4 调用流程、编译运行说明。

## 核心实现位置

- `../sm4_channel.h`: 对外接口和结果结构体。
- `../sm4_channel.cpp`: SM2 KEX 封装、SM3 派生、SM4-CBC 加解密实现。

## 编译运行

从 `src` 目录执行：

```sh
make secure_message_demo
./secure_message_demo
```

或者：

```sh
make run-secure-message-demo
```

菜单中的正确性测试和性能评测会直接调用当前目录的 Makefile：

```sh
make run-pornin-inv-asm USE_FULL_INV_ASM=1
make perf-pornin-asm USE_FULL_INV_ASM=1
```

因此建议始终在 `src` 目录启动 `./secure_message_demo`。

## 分层说明

本模块只调用已有 SM2 密钥交换上层接口。当前 `secure_message_demo` 构建目标链接的是 Montgomery + fp ASM + `src/asm/pornin_full_inv.S` 的优化路径，不修改这些底层优化逻辑。
