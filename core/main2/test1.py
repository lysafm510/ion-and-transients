from math import log, pi

import numpy as np
import pandas as pd

from core.main.constant import F, Z, T, R, DT, C

# coeff = - R * T / (Z * F)
# # fn = 0.9211692111937463
# fn = 0.48224551805931404
# CCACYTO = 0.0001
# nernst = coeff * log(fn / CCACYTO)
#
# print(nernst)
#
# # i_ca_ryr = 2.5319073718753766
# i_ca_ryr = 0.8949798651605706
# # b = (i_ca_ryr * 10 ** (-12)) * DT / (C * 10 ** (-20))
#
# b = (i_ca_ryr * DT) / (C * 10 ** (-6))
#
# print(b)
# #
# c = 1 - b/nernst
# print(c)


# i_ca_ryr = 2.5319073718753766
# vm = 0
# # 算出电流之后，继续计算当前膜电位
# area = pi * 300 ** 2
# vm = (vm - i_ca_ryr * DT / (C * area)) * (10.0 ** 11)
# print(vm)
# # 计算平衡电位
# c_ca_store = 0.9211692111937463
# CCAMYO = 0.0001
# nernst = - R * T / (Z * F) * log(c_ca_store / CCAMYO) * (10 ** 3)
# print(nernst)
# print(vm / nernst)
# # 计算DCARYR系数
# f = (1 - vm / nernst)
# print(f)

RELEASE_TIMES = 2 * 10 ** -2
RELEASE_STEP = int((RELEASE_TIMES + (10 ** -6)) / DT)
print(RELEASE_STEP)