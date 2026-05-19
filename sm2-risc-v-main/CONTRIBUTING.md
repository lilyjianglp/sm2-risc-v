# 贡献指南

感谢你关注本项目。本项目面向 RISC-V 平台上的国密密码工程优化，核心是 **SM2 密钥交换及其底层有限域运算优化**，并在上层提供 **SM4 安全通信 CLI 演示** 和 **前端展示页面**。

本指南说明如何提交 issue、创建分支、提交 PR、运行测试，以及贡献代码时需要遵守的基本规范。

## 项目定位

本项目主要包含三类内容：

```text
1. SM2 RISC-V 优化后端
2. SM4 安全通信 CLI 演示
3. React 前端展示页面
```

推荐保持如下分层：

```text
CLI / frontend display layer
-> SM2 key exchange protocol layer
-> SM2 scalar multiplication layer
-> Montgomery finite-field backend
-> RISC-V assembly optimization
-> SM3/KDF key derivation
-> SM4-CBC secure message demo
```

新增功能应优先调用已有上层接口，不建议绕过协议层直接改动底层汇编细节。

## 目录说明

典型目录结构如下：

```text
.
├── README.md
├── CONTRIBUTING.md
├── frontend/
└── sm2-risc-v-main/
    └── src/
        ├── asm/
        ├── demo/
        ├── Makefile
        ├── sm2_kex*.c
        ├── sm2_scalar*.c
        ├── fp*.c
        └── sm4_channel.*
```

主要目录说明：

| 路径 | 作用 |
|---|---|
| `sm2-risc-v-main/src/` | RISC-V 后端核心代码和 Makefile |
| `sm2-risc-v-main/src/asm/` | RISC-V 汇编优化代码 |
| `sm2-risc-v-main/src/demo/` | CLI 交互演示入口 |
| `frontend/` | React 前端展示页面 |
| `frontend/src/demoData.js` | 前端展示用的板端样例数据 |

## 提交 Issue

如果你发现 bug、编译失败、测试失败、性能异常，或希望提出功能建议，请提交 issue。

Issue 建议包含：

```text
1. 问题描述
2. 复现步骤
3. 运行环境
4. 执行命令
5. 实际输出
6. 期望结果
```

运行环境建议写清楚：

```text
Board:
SoC:
OS:
Compiler:
Branch:
Build command:
Run command:
```

例如：

```text
Board: Muse Pi Pro
SoC: SpacemiT K1
OS: Linux RISC-V 64-bit
Compiler: gcc/g++
Command: make secure_message_demo
```

## Fork 与分支

建议使用 Fork + Branch + Pull Request 的方式参与开发。

```bash
git clone https://github.com/<your-name>/sm2-risc-v.git
cd sm2-risc-v
git remote add upstream https://github.com/lilyjianglp/sm2-risc-v.git
git remote -v
```

不要直接在 `main` 分支上开发。请为每个任务创建独立分支：

```bash
git checkout -b feature/sm4-secure-demo
```

分支命名建议：

| 类型 | 示例 |
|---|---|
| 新功能 | `feature/sm4-secure-demo` |
| bug 修复 | `fix/build-error` |
| 文档 | `docs/update-readme` |
| 测试 | `test/sm2-kex-correctness` |
| 性能 | `bench/montgomery-backend` |
| 构建 | `build/makefile-cleanup` |

## 提交前检查

提交 PR 前，请根据修改范围运行对应检查。

### CLI / 后端演示

如果修改了 SM4 安全通信 demo、CLI 菜单或 `sm4_channel.*`，请运行：

```bash
cd sm2-risc-v-main/src
make secure_message_demo
./secure_message_demo
```

重点检查：

```text
[SM2 KEX] Shared key generated successfully.
[SM2 KEX] KA == KB: true
[PASS] Secure message demo success.
```

### 汇编 / Montgomery / SM2 后端

如果修改了 SM2、Montgomery、有限域运算或 RISC-V 汇编相关代码，请尽量运行正确性测试：

```bash
cd sm2-risc-v-main/src
make run-pornin-inv-asm USE_FULL_INV_ASM=1
```

如果修改了性能关键路径，请运行 benchmark：

```bash
cd sm2-risc-v-main/src
make perf-pornin-asm USE_FULL_INV_ASM=1
```

如需开启用户态 perf 访问，可在板端执行：

```bash
sudo sysctl -w kernel.perf_user_access=2
```

### 前端

如果修改了 `frontend/`，请运行：

```bash
cd frontend
npm install
npm run build
```

本地预览：

```bash
npm run dev -- --host 127.0.0.1 --port 3000
```

Windows PowerShell 可使用：

```powershell
npm.cmd install
npm.cmd run build
npm.cmd run dev -- --host 127.0.0.1 --port 3000
```

## Pull Request 要求

提交 PR 时，请确认：

- 分支基于最新主分支；
- PR 标题能概括修改内容；
- PR 描述写清楚修改原因；
- 已列出测试命令和测试结果；
- 没有提交本地临时文件、编译产物或个人绝对路径；
- Makefile 中使用相对路径，不写死个人机器路径。

PR 描述建议格式：

```markdown
## Summary

说明本次修改内容。

## Changes

- 修改 1
- 修改 2
- 修改 3

## Test

```bash
cd sm2-risc-v-main/src
make secure_message_demo
./secure_message_demo
```

## Notes

说明是否影响 SM2、Montgomery、RISC-V ASM、Makefile 或前端展示。
```

## 代码风格

### C / C++

- 保持现有命名和文件组织风格；
- 新增功能优先复用已有接口；
- 不要在核心密码路径中加入无关调试输出；
- 不要大规模重构与当前任务无关的代码；
- 安全通信流程应保持清晰分层；
- 修改 SM2、Montgomery 或汇编相关代码后必须补充测试说明。

### RISC-V 汇编

- 不要随意修改底层汇编优化代码；
- 修改汇编前应明确对应的 C 参考实现；
- 修改后需要运行正确性测试；
- 性能优化 PR 建议附带 benchmark 数据；
- 避免提交只有格式变化、没有实际逻辑意义的汇编改动。

### Makefile

- 不要写死个人本地路径，例如：

```text
C:\Users\...
/home/<your-name>/...
```

- 优先使用相对路径，例如：

```text
asm/pornin_full_inv.S
demo/secure_message_demo.cpp
sm4_channel.cpp
```

- 编译器应允许外部覆盖，例如：

```bash
make CC=riscv64-linux-gnu-gcc CXX=riscv64-linux-gnu-g++
```

### 前端

- 前端只作为展示层；
- 不在前端重新实现真实密码算法；
- 板端真实样例数据放在 `frontend/src/demoData.js`；
- 展示内容应与仓库实际代码保持一致；
- 如果某个优化实现尚未开源，不要在页面中写成“已接入”。

### 文档

- 命令应能直接复制运行；
- 路径应与仓库结构一致；
- 已完成内容和计划内容要分开写；
- benchmark 数据应注明测试平台和测试命令；
- 不写无法验证的性能结论。

## Commit Message 规范

推荐格式：

```text
<type>: <description>
```

常用类型：

| 类型 | 含义 |
|---|---|
| `feat` | 新增功能 |
| `fix` | 修复问题 |
| `docs` | 文档修改 |
| `test` | 测试相关 |
| `bench` | 性能评测相关 |
| `perf` | 性能优化 |
| `build` | 构建系统修改 |
| `refactor` | 重构 |
| `chore` | 维护类修改 |

示例：

```bash
git commit -m "feat: add SM4 secure message demo"
git commit -m "docs: update SM2 backend usage guide"
git commit -m "bench: add RISC-V board benchmark result"
git commit -m "build: use relative paths in Makefile"
```

## 不建议提交的内容

请不要提交：

- `node_modules/`
- `frontend/dist/`
- 本地编译生成的可执行文件；
- `.dump`、`perf.data`、`perf.data.old` 等临时文件；
- IDE 缓存；
- 个人绝对路径；
- 未确认来源的 benchmark 数据；
- 与当前 PR 无关的大规模格式化修改。

如确实需要保留某些生成文件，请在 PR 中说明原因。

## 贡献方向

欢迎贡献以下内容：

```text
1. RISC-V 平台适配
2. SM2 密钥交换优化
3. Montgomery 有限域运算优化
4. RISC-V 汇编正确性测试
5. benchmark 数据与分析
6. SM4 安全通信 demo 改进
7. 前端展示优化
8. README / docs 文档完善
```

如果不确定某项修改是否适合提交，请先创建 issue 讨论。
