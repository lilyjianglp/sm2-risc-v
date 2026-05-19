export const boardDemoRuns = [
  {
    message: 'hello im Aria',
    sharedKey: '310bb44bb38989e769309604c44bee80b4c8cc9d6fdde84af49d7d4712857eb0',
    sm4Key: '3b5342a6ab9c60d09486788aae0192ba',
    iv: '311ffaf89730283aeb052da2517a8ffa',
    ciphertext: 'cf532250028e91e6dbe4b990d6f43c0b',
    decryptedText: 'hello im Aria',
    kexSuccess: true,
    decryptSuccess: true,
    source: 'RISC-V board CLI run'
  },
  {
    message: 'risc-v sm2 kex ready',
    sharedKey: '310bb44bb38989e769309604c44bee80b4c8cc9d6fdde84af49d7d4712857eb0',
    sm4Key: '3b5342a6ab9c60d09486788aae0192ba',
    iv: '0dded5e2efad25fb5bbd326bc01b6f3f',
    ciphertext: '926a07a548968a08376cd6d0459ca3b8162b412bc7db7a69fca3e5f630f098c3',
    decryptedText: 'risc-v sm2 kex ready',
    kexSuccess: true,
    decryptSuccess: true,
    source: 'RISC-V board CLI run'
  },
  {
    message: 'secure channel message',
    sharedKey: '310bb44bb38989e769309604c44bee80b4c8cc9d6fdde84af49d7d4712857eb0',
    sm4Key: '3b5342a6ab9c60d09486788aae0192ba',
    iv: '399ebc222e1609fbdd12473c122cba72',
    ciphertext: '1a8161c90f8571e03b28f2da40fc456f5f32377a35b2a9ab0b25d9b6318a1797',
    decryptedText: 'secure channel message',
    kexSuccess: true,
    decryptSuccess: true,
    source: 'RISC-V board CLI run'
  },
  {
    message: 'SM4 demo on Muse Pi Pro',
    sharedKey: '310bb44bb38989e769309604c44bee80b4c8cc9d6fdde84af49d7d4712857eb0',
    sm4Key: '3b5342a6ab9c60d09486788aae0192ba',
    iv: '2c9de4b02e7ab87f9259af2af35286e0',
    ciphertext: '0c5fd99b9a4407f442aaa6d3bf86ee1819b9485c90681b4a69f4ad21d6dd79c5',
    decryptedText: 'SM4 demo on Muse Pi Pro',
    kexSuccess: true,
    decryptSuccess: true,
    source: 'RISC-V board CLI run'
  }
];

export const sampleMessages = [
  'hello im Aria',
  'risc-v sm2 kex ready',
  'secure channel message',
  'SM4 demo on Muse Pi Pro'
];
