from math import log
from constant import DCARYR, DT, C, R, T, Z, Far, K

ca_in1 = 0.719930926
ca_in2 = 0.63471069
ca_out = 0.0001
nernst1 = - R * T / (Z * Far) * log(ca_in1 / ca_out) * (10.0 ** 3)
nernst2 = - R * T / (Z * Far) * log(ca_in2 / ca_out) * (10.0 ** 3)
print("平均浓度", nernst1)
print("出流口浓度", nernst2)
