    .section .text
    .global read_cycle

read_cycle:
    rdcycle a0   # 使用 rdcycle 获取时钟周期数，结果保存在 a0 寄存器
    ret          # 返回
