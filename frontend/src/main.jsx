import React, { useMemo, useState } from 'react';
import { createRoot } from 'react-dom/client';
import { boardDemoRuns, sampleMessages } from './demoData.js';
import './styles.css';

const pages = [
  'SM2 优化概览',
  'SM2 密钥交换',
  'SM4 安全信道',
  'RISC-V 后端链路',
  '性能基准评测',
  '运行指南'
];

const fallbackRun = boardDemoRuns[0];

const coreMetrics = [
  { name: 'Montgomery 模乘', value: '330', unit: 'cycles', category: 'SM2 底层', tone: 'blue' },
  { name: 'Montgomery 模平方', value: '298', unit: 'cycles', category: 'SM2 底层', tone: 'blue' },
  { name: 'Pornin 全模逆 ASM', value: '17,907', unit: 'cycles', category: 'SM2 底层', tone: 'green' },
  { name: 'SM2 KEX 发起方', value: '2,846,039', unit: 'cycles', category: 'SM2 密钥交换', tone: 'amber' },
  { name: 'SM2 KEX 响应方', value: '3,128,943', unit: 'cycles', category: 'SM2 密钥交换', tone: 'amber' },
  { name: 'IPC 指令吞吐率', value: '1.36', unit: 'insn/cycle', category: 'CPU 特性', tone: 'green' }
];

const sm2Optimizations = [
  { title: 'Montgomery 有限域后端', desc: '承接 SM2 曲线运算中的模乘、模平方与模约减。', highlight: '核心优化' },
  { title: 'Pornin 全模逆汇编', desc: '将 fp_mont_inv 关键路径接入 RISC-V 汇编实现。', highlight: '核心优化' },
  { title: 'SM2 标量乘法', desc: '服务密钥交换中的点乘与验证流程。', highlight: '关键路径' },
  { title: 'SM2 密钥交换', desc: '完成发起方与响应方共享密钥协商，并验证 KA == KB。', highlight: '协议层' }
];

const backendPipeline = [
  ['CLI / 前端展示层', '用户交互与结果展示'],
  ['SM2 密钥交换协议层', '协商共享密钥 KA / KB，验证一致性'],
  ['SM2 标量乘法层', '执行点乘、倍点与验证相关运算'],
  ['Montgomery 有限域后端', '模乘 / 模平方 / 模约减优化'],
  ['Pornin 全模逆汇编', 'fp_mont_inv ASM 关键路径加速'],
  ['SM3 密钥派生层', '从 SM2 共享密钥派生 SM4 会话密钥'],
  ['SM4-CBC 安全信道', '消息加密与解密验证']
];

function hexFromMessage(message) {
  let x = 0x811c9dc5;
  for (let i = 0; i < message.length; i += 1) {
    x ^= message.charCodeAt(i);
    x = Math.imul(x, 0x01000193) >>> 0;
  }
  const words = [];
  for (let i = 0; i < 4; i += 1) {
    x ^= x << 13;
    x ^= x >>> 17;
    x ^= x << 5;
    words.push((x >>> 0).toString(16).padStart(8, '0'));
  }
  return words.join('');
}

function maskHex(value, keep = 8) {
  const text = String(value ?? '');
  if (text.length <= keep * 2 + 3) {
    return text;
  }
  return `${text.slice(0, keep)}...${text.slice(-keep)}`;
}

function sourceText(source) {
  if (!source) return '板端验证数据';
  if (source.includes('frontend preview')) return '前端预览数据';
  if (source.includes('board')) return 'RISC-V 板端实测数据';
  return source;
}

function App() {
  const [message, setMessage] = useState('hello im Aria');
  const [pageIndex, setPageIndex] = useState(0);
  const activePage = pages[pageIndex];
  const previewCiphertext = useMemo(() => hexFromMessage(message), [message]);
  const storedRun = useMemo(
    () => boardDemoRuns.find((run) => run.message === message.trim()),
    [message]
  );
  const displayRun = storedRun ?? {
    ...fallbackRun,
    ciphertext: previewCiphertext,
    decryptedText: message,
    source: 'frontend preview only'
  };

  const prevPage = () => setPageIndex((index) => Math.max(index - 1, 0));
  const nextPage = () => setPageIndex((index) => Math.min(index + 1, pages.length - 1));

  return (
    <main className="app-shell">
      <aside className="sidebar">
        <div className="brand-lockup">
          <div className="brand-mark">SM2</div>
          <div>
            <p className="eyebrow">RV-GM-Secure</p>
            <h1>
              国密安全通信
              <span>SM2 优化 + SM4 安全信道</span>
            </h1>
          </div>
        </div>

        <nav className="side-nav" aria-label="页面导航">
          {pages.map((page, index) => (
            <button
              className={index === pageIndex ? 'active' : ''}
              key={page}
              type="button"
              onClick={() => setPageIndex(index)}
            >
              {String(index + 1).padStart(2, '0')} {page}
            </button>
          ))}
        </nav>

        <div className="board-badge">
          <span className="status-dot" />
          <div>
            <strong>RISC-V Muse Pi Pro 板端已验证</strong>
            <small>SpacemiT K1 / SM2 KEX 优化后端</small>
          </div>
        </div>
      </aside>

      <section className="workspace">
        <header className="topbar">
          <div>
            <p className="eyebrow">
              {activePage === 'SM2 密钥交换' || activePage === 'SM4 安全信道'
                ? 'SM2 KEX → SM3 KDF → SM4-CBC 完整链路'
                : 'RISC-V 平台 · 汇编优化后端'}
            </p>
            <h2>{activePage}</h2>
          </div>
          <div className="pager">
            <button type="button" onClick={prevPage} disabled={pageIndex === 0}>
              上一页
            </button>
            <span>{pageIndex + 1} / {pages.length}</span>
            <button type="button" onClick={nextPage} disabled={pageIndex === pages.length - 1}>
              下一页
            </button>
          </div>
        </header>

        {activePage === 'SM2 优化概览' && (
          <section className="panel overview-panel">
            <div className="section-heading">
              <div>
                <p className="eyebrow">项目核心</p>
                <h3>RISC-V 平台 SM2 密钥交换底层优化 + SM4 安全通信演示</h3>
              </div>
              <span className="panel-tag">优化已验证</span>
            </div>

            <div className="optimization-highlights">
              <div className="optimization-header">
                <span className="badge">核心优化</span>
                <strong>SM2 在 RISC-V 平台上的关键路径优化</strong>
              </div>
              <div className="optimization-grid">
                {sm2Optimizations.map((opt) => (
                  <article key={opt.title} className="optimization-card">
                    <div className="opt-tag">{opt.highlight}</div>
                    <strong>{opt.title}</strong>
                    <span>{opt.desc}</span>
                  </article>
                ))}
              </div>
            </div>

            <div className="overview-grid">
              <article>
                <strong>SM2 密钥交换优化链路</strong>
                <span>
                  SM2 KEX → sm2_kex_mont.c → sm2_scalar_mont.c → Montgomery 后端 → Pornin 全模逆 ASM。
                </span>
              </article>
              <article>
                <strong>完整安全通信</strong>
                <span>
                  SM2 共享密钥 → SM3 KDF → SM4 会话密钥 → SM4-CBC 加解密。
                </span>
              </article>
              <article>
                <strong>汇编优化关键路径</strong>
                <span>
                  fp_mont_mul / fp_mont_sqr / fp_mont_inv 连接 RISC-V 汇编优化实现。
                </span>
              </article>
              <article>
                <strong>板端实测性能</strong>
                <span>
                  SM2 KEX 发起方约 2.85M cycles，响应方约 3.13M cycles，IPC 1.36。
                </span>
              </article>
            </div>

            <div className="hero-readout">
              <span>验证状态</span>
              <strong>[PASS] SM2 优化后端 + SM4 安全通信完整链路验证通过</strong>
              <div className="readout-stats">
                <span>Montgomery 模乘: 330 cycles</span>
                <span>Pornin 模逆 ASM: 17,907 cycles</span>
                <span>KA == KB: true</span>
              </div>
            </div>
          </section>
        )}

        {activePage === 'SM2 密钥交换' && (
          <section className="panel sm2-panel">
            <div className="section-heading">
              <div>
                <p className="eyebrow">SM2 优化密钥交换</p>
                <h3>共享密钥协商 · 一致性验证</h3>
              </div>
              <span className="panel-tag">
                {storedRun ? '板端实测数据' : '前端预览模式'}
              </span>
            </div>

            <div className="sm2-kex-grid">
              <div className="sm2-optimization-note">
                <div className="note-header">
                  <span>SM2 优化后端</span>
                </div>
                <div className="note-content">
                  <code>Montgomery 有限域后端 + Pornin 全模逆 RISC-V 汇编</code>
                  <div className="opt-badges">
                    <span>fp_mont_mul: 330 cycles</span>
                    <span>fp_mont_inv ASM: 17,907 cycles</span>
                  </div>
                </div>
              </div>

              <div className="message-console">
                <label htmlFor="message">测试消息：后续使用 SM2 共享密钥派生出的 SM4 密钥加密</label>
                <textarea
                  id="message"
                  value={message}
                  onChange={(event) => setMessage(event.target.value)}
                  spellCheck="false"
                />

                <div className="sample-row">
                  <div className="sample-heading">
                    <strong>板端真实测试样例</strong>
                    <span>点击样例展示 RISC-V 板端实测的 SM2 KEX + SM4 加密数据。</span>
                  </div>
                  {sampleMessages.map((sample) => (
                    <button key={sample} type="button" onClick={() => setMessage(sample)}>
                      {sample}
                    </button>
                  ))}
                </div>
              </div>

              <div className="sm2-result-readout">
                <div className="result-header">
                  <strong>SM2 密钥交换结果</strong>
                  <span className="success-badge">优化后端</span>
                </div>
                <Readout label="SM2 共享密钥摘要" value={maskHex(displayRun.sharedKey)} />
                <Readout
                  label="KA == KB 一致性验证"
                  value={storedRun ? 'true，板端验证通过' : 'true，预览模式'}
                  plain
                  highlight
                />
                <Readout
                  label="SM2 KEX 状态"
                  value={storedRun ? '完成，RISC-V 板端实测' : '前端预览'}
                  plain
                />
              </div>
            </div>

            <div className="sm2-to-sm4-note">
              <span>SM2 共享密钥将通过 SM3 KDF 派生为 128-bit SM4 会话密钥，用于安全信道。</span>
            </div>
          </section>
        )}

        {activePage === 'SM4 安全信道' && (
          <section className="panel channel-panel">
            <div className="section-heading">
              <div>
                <p className="eyebrow">SM2 共享密钥 → SM3 KDF → SM4-CBC</p>
                <h3>基于 SM2 密钥交换的安全通信演示</h3>
              </div>
              <span className="panel-tag">
                {storedRun ? '板端实测数据' : '前端预览模式'}
              </span>
            </div>

            <div className="channel-grid">
              <div className="crypto-readout">
                <Readout label="数据来源" value={sourceText(displayRun.source)} plain />
                <div className="readout-divider">
                  <span>源自 SM2 密钥交换</span>
                </div>
                <Readout label="SM2 共享密钥摘要" value={maskHex(displayRun.sharedKey)} />
                <Readout label="SM3 KDF 派生的 SM4 会话密钥" value={maskHex(displayRun.sm4Key)} />
                <Readout label="SM4-CBC 初始化向量 IV" value={displayRun.iv} />
                <Readout label="密文 Ciphertext" value={displayRun.ciphertext} />
                <Readout label="解密后明文" value={displayRun.decryptedText || ' '} plain />
              </div>

              <div className="kdf-flow">
                <div className="flow-header">
                  <strong>密钥派生流程</strong>
                </div>
                <div className="flow-steps">
                  <FlowStep index="1" title="SM2 共享密钥" caption="KA / KB 256-bit" />
                  <div className="flow-arrow">↓</div>
                  <FlowStep index="2" title="SM3 KDF" caption="密钥派生函数" />
                  <div className="flow-arrow">↓</div>
                  <FlowStep index="3" title="SM4 会话密钥" caption="128-bit" />
                  <div className="flow-arrow">↓</div>
                  <FlowStep index="4" title="SM4-CBC" caption="加密 / 解密消息" />
                </div>
              </div>
            </div>

            <CliTranscript run={displayRun} inputMessage={message} stored={Boolean(storedRun)} />
          </section>
        )}

        {activePage === 'RISC-V 后端链路' && (
          <section className="panel">
            <div className="section-heading">
              <div>
                <p className="eyebrow">汇编优化调用链</p>
                <h3>SM2 KEX → Montgomery → Pornin ASM 完整后端路径</h3>
              </div>
              <span className="panel-tag">RISC-V ASM</span>
            </div>

            <div className="pipeline">
              {backendPipeline.map(([title, caption], index) => {
                const highlighted = title.includes('Pornin') || title.includes('Montgomery');
                return (
                  <React.Fragment key={title}>
                    <div className={`pipeline-node ${highlighted ? 'highlight-node' : ''}`}>
                      <span>{String(index + 1).padStart(2, '0')}</span>
                      <strong>{title}</strong>
                      <small>{caption}</small>
                      {highlighted && <div className="asm-badge">RISC-V ASM</div>}
                    </div>
                    {index < backendPipeline.length - 1 && <div className="pipeline-edge" />}
                  </React.Fragment>
                );
              })}
            </div>

            <div className="asm-files">
              <strong>汇编优化源文件</strong>
              <div className="asm-list">
                <code>asm/pornin_full_inv.S</code>
                <code>asm/fp_mont_mul.S</code>
                <code>asm/fp_mont_sqr.S</code>
                <code>asm/fn_mont.S</code>
              </div>
            </div>
          </section>
        )}

        {activePage === '性能基准评测' && (
          <section className="panel">
            <div className="section-heading">
              <div>
                <p className="eyebrow">RISC-V 板端实测数据</p>
                <h3>SpacemiT K1 / Muse Pi Pro 性能基准</h3>
              </div>
              <span className="panel-tag">perf stat + rdcycle</span>
            </div>

            <MetricSection
              title="SM2 底层运算优化成果"
              metrics={coreMetrics.filter((metric) => metric.category === 'SM2 底层')}
            />
            <MetricSection
              title="SM2 密钥交换性能"
              metrics={coreMetrics.filter((metric) => metric.category === 'SM2 密钥交换')}
            />
            <MetricSection
              title="系统综合指标"
              metrics={coreMetrics.filter((metric) => metric.category === 'CPU 特性')}
            />

            <div className="perf-strip">
              <span>用户态周期数: 9,946,462,245</span>
              <span>用户态指令数: 13,498,184,886</span>
              <span>总耗时: 6.222s</span>
              <span className="highlight">SM2 KEX 优化链路板端实测通过</span>
            </div>
          </section>
        )}

        {activePage === '运行指南' && (
          <section className="panel command-panel">
            <div className="section-heading">
              <div>
                <p className="eyebrow">板端运行命令</p>
                <h3>编译、运行与性能测试入口</h3>
              </div>
            </div>

            <div className="command-grid">
              <div className="command-block">
                <div className="command-title">编译与运行 SM2 + SM4 演示</div>
                <CodeBlock lines={['cd src', 'make secure_message_demo', './secure_message_demo']} />
              </div>
              <div className="command-block">
                <div className="command-title">测试 Pornin 模逆汇编优化</div>
                <CodeBlock
                  lines={[
                    'make run-pornin-inv-asm USE_FULL_INV_ASM=1',
                    'make perf-pornin-asm USE_FULL_INV_ASM=1'
                  ]}
                />
              </div>
              <div className="command-block full-width">
                <div className="command-title">CLI 菜单入口</div>
                <CodeBlock
                  lines={[
                    '1. Run SM2 key exchange demo',
                    '2. Run SM4 secure message demo',
                    '3. Run correctness test',
                    '4. Run performance benchmark'
                  ]}
                />
              </div>
            </div>
          </section>
        )}
      </section>
    </main>
  );
}

function MetricSection({ title, metrics }) {
  return (
    <div className="metric-section">
      <div className="metric-section-title">
        <span>{title}</span>
      </div>
      <div className="metric-grid">
        {metrics.map((metric) => (
          <article className={`metric-card ${metric.tone}`} key={metric.name}>
            <span>{metric.name}</span>
            <strong>{metric.value}</strong>
            <small>{metric.unit}</small>
          </article>
        ))}
      </div>
    </div>
  );
}

function FlowStep({ index, title, caption }) {
  return (
    <div className="flow-step">
      <span>{index}</span>
      <code>{title}</code>
      <small>{caption}</small>
    </div>
  );
}

function Readout({ label, value, plain = false, highlight = false }) {
  return (
    <div className={`readout ${highlight ? 'highlight-readout' : ''}`}>
      <span>{label}</span>
      <code className={plain ? 'plain-text' : ''}>{value}</code>
    </div>
  );
}

function CliTranscript({ run, inputMessage, stored }) {
  const safeMessage = inputMessage || '<空消息>';
  const lines = [
    '=========================================',
    ' RV-GM-Secure: RISC-V SM2 + SM4 Backend',
    '=========================================',
    '2. Run SM4 secure message demo (based on SM2 KEX)',
    'Please select: 2',
    '',
    'Please input message:',
    safeMessage,
    '',
    '[SM2 KEX] Running optimized SM2 key exchange...',
    '[SM2 KEX] Montgomery backend + Pornin full inversion ASM',
    '[SM2 KEX] Shared key generated successfully.',
    `[SM2 KEX] KA == KB: true${stored ? ' (RISC-V board verified)' : ' (preview)'}`,
    '',
    '[KDF] Deriving SM4 session key from SM2 shared key (SM3 KDF)...',
    `[KDF] SM4 key digest: ${maskHex(run.sm4Key)}`,
    `[SM4-CBC] IV: ${run.iv}`,
    '[SM4-CBC] Encrypting message...',
    `[SM4-CBC] Ciphertext: ${run.ciphertext}`,
    '',
    '[SM4-CBC] Decrypting message...',
    `[SM4-CBC] Plaintext: ${run.decryptedText || safeMessage}`,
    '',
    stored
      ? '[PASS] Secure message demo success on RISC-V board.'
      : '[PREVIEW] 当前输入未匹配板端实测记录，展示前端预览结果。'
  ];

  return (
    <div className="cli-transcript">
      <div className="terminal-bar">
        <span />
        <strong>RISC-V 板端 CLI 运行日志预览</strong>
        <em>SM2 optimized backend</em>
      </div>
      <pre>
        {lines.map((line, index) => (
          <code key={`${index}-${line}`}>{line}</code>
        ))}
      </pre>
    </div>
  );
}

function CodeBlock({ lines }) {
  return (
    <pre className="code-block">
      {lines.map((line) => (
        <code key={line}>{line}</code>
      ))}
    </pre>
  );
}

createRoot(document.getElementById('root')).render(<App />);
