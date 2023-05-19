from math import log, sqrt

import numpy as np

from constant import B_OUTFLOW, F, Z, T, R
from grid_info import NP, NE, NOD, NPOCH, PX, PY


def cal_nernst_potential(fn, CCACYTO):
    nernst = np.zeros(NP, float)
    count = 0
    for i in range(0, NP):
        if NPOCH[i] == B_OUTFLOW:
            count = count + 1
            coeff = - R * T / (Z * F)
            nernst[i] = coeff * log(fn[i] / CCACYTO)
    avg_nernst = nernst.sum() / count

    return avg_nernst
