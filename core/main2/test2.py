import math

from core.main.constant import DT

Vm = 0
Cm = 2.0
R = 8.314
T = 310.0
U = 1.610217733 * 10 ** -19  # 一个电荷的常数
M = 6.0221367 * 10 ** 23  # 摩尔常数，单位
cca_out = 0.0001  # 初始膜外钙离子浓度 有用到初始膜内值？

K_RYR = 6.5 * 10 ** 7  # *0.25？

ICa = 2.5319073718753766
Vm = (ICa * DT / Cm + Vm) * (10.0 ** -3)


# 平衡状态下
cca_in = 0.9211692111937463  # 边界膜内钙离子浓度
c = math.log(cca_in / cca_out, math.e)
Vj = R * T / (2 * U * M) * c * (10.0 ** 3)  # 还加负号吗，Vm值是正的吧，虽然膜内电压是负的

f = Vm / Vj

print("vm, ", Vm)
print("vj, ", Vj)
print("f, ", f)
print("(1-f), ", (1 - f))
DCARYR = (1 - f) * K_RYR
