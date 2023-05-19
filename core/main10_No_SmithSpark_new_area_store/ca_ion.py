from math import log, pi
from constant import DCARYR, DT, C, R, T, Z, Far, K, R_JSR, H_JSR
from blink_coeff_matrix import TOTAL_AREA


def get_d_ca_ryr(vm, i_ca_ryr, i_ca_fsr, ca_in, ca_out):
    # 计算膜电位
    ica = i_ca_ryr - i_ca_fsr
    # 计算终池表面积
    area = TOTAL_AREA * 2 + 2 * pi * R_JSR * H_JSR
    vm_new = vm - (1 - K) * ica * DT / (C * area) * (10.0 ** 11)

    # 计算能斯特电位
    nernst = - R * T / (Z * Far) * log(ca_in / ca_out) * (10.0 ** 3)

    d_ca_ryr = (1 - vm_new / nernst) * DCARYR

    return vm_new, nernst, d_ca_ryr
