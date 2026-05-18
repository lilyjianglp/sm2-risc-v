# RV-GM-Secure Frontend

这是一个用于展示“国密安全通信后端 / RISC-V 密码工程优化”的 React + Vite 前端页面。

页面重点展示：

- SM2 密钥交换协议层：通过上层 SM2 KEX 接口生成 KA / KB，并验证双方共享密钥一致；
- Montgomery 有限域后端：支撑 SM2 曲线运算中的模乘、模平方和模约减；
- Pornin 全模逆汇编优化：用于优化 `fp_mont_inv` 关键路径；
- SM3/KDF 密钥派生：从 SM2 shared key 派生 128-bit SM4 会话密钥；
- SM4-CBC 安全信道：对用户消息进行加密、解密和一致性验证；
- RISC-V 板端性能数据：展示 Muse Pi Pro 上的周期数、指令数和 IPC 结果。

## 运行方式

```sh
cd frontend
npm install
npm run dev
```

如果 Windows PowerShell 禁止执行 `npm.ps1`，使用：

```powershell
npm.cmd install
npm.cmd run dev -- --host 127.0.0.1 --port 3000
```

浏览器打开：

```text
http://127.0.0.1:3000/
```

## 数据来源

当前前端是展示型页面，不直接调用 RISC-V 板端程序。

真实后端入口仍然是：

```sh
cd src
make secure_message_demo
./secure_message_demo
```

页面中的板端样例来自已经在 RISC-V Muse Pi Pro 上跑通的 CLI 输出。输入框支持自由输入：

- 如果输入内容与板端样例完全一致，页面展示真实板端记录；
- 如果输入新内容，页面进入前端预览模式，只用于展示交互效果；
- 真正的 SM2 KEX、SM3/KDF、SM4-CBC 加解密仍以 `src/secure_message_demo` 的运行结果为准。

## 添加板端样例

在板子上多跑几组 CLI 后，把结果追加到：

```text
frontend/src/demoData.js
```

字段格式：

```js
{
  message: 'your input',
  sharedKey: '...',
  sm4Key: '...',
  iv: '...',
  ciphertext: '...',
  decryptedText: 'your input',
  kexSuccess: true,
  decryptSuccess: true,
  source: 'RISC-V board CLI run'
}
```
