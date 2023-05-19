from math import log
from constant import DCARYR, DT, C, R, T, Z, Far
from blink_coeff_matrix import TOTAL_AREA


def get_d_ca_ryr(vm, ica, ca_in, ca_out):
    # vm2 = vm - ica * DT / (C * TOTAL_AREA) * (10.0 ** 11)  # 计算膜电位
    vm2 = vm - ica * DT / C * (10.0 ** -3)
    nernst = - R * T / (Z * Far) * log(ca_in / ca_out) * (10.0 ** 3)  # 计算能斯特电位
    f = 1 - vm2 / nernst
    d_ca_ryr = f * DCARYR
    print(vm2)
    print(nernst)
    print(f)
    print(d_ca_ryr)
    return vm2, d_ca_ryr
