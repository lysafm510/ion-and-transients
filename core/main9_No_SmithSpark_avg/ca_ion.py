from math import log
from constant import DCARYR, DT, C, R, T, Z, Far, K
from blink_coeff_matrix import TOTAL_AREA


def get_d_ca_ryr(vm, i_ca_ryr, i_ca_fsr, ca_in, ca_out):

    # 计算膜电位
    ica = i_ca_ryr - i_ca_fsr
    print()
    print("ica", ica)

    vm_new = vm - (1 - K) * ica * DT / (C * TOTAL_AREA) * (10.0 ** 11)

    print("ca_in", ca_in)
    print("ca_out", ca_out)
    # 计算能斯特电位
    nernst = - R * T / (Z * Far) * log(ca_in / ca_out) * (10.0 ** 3)

    f = 1 - vm_new / nernst
    d_ca_ryr = f * DCARYR

    print("vm: ", vm_new)
    print("nernst: ", nernst)
    print("DCARYR: ", d_ca_ryr)
    print()
    return vm_new, nernst, d_ca_ryr
