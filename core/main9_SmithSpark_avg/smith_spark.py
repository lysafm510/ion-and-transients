from math import pi
import numpy as np
from constant import NR, KPCAF, BCAF, KMCAF, KPCAM, BCAM, KMCAM, KMTRC, KMSRM, KMSLM, BTRC, BSRM, BSLM, KPSLM, \
    KPSRM, KPTRC, DCACYT, DT, DR, DCAF
from grid_info import RADIUS


def cal_dye(j_dye, c_ca_cyt, c_caf):
    # 计算ca扩散方程的dye项,c_ca_cyt，c_caf用n时刻的值

    for i in range(0, NR):
        j_dye[i] = -KPCAF * c_ca_cyt[i] * (BCAF - c_caf[i]) + KMCAF * c_caf[i]
    return j_dye


def cal_buffers(j_buffers, c_ca_cyt, c_cam, c_trc, c_srm, c_slm):
    # 计算ca扩散方程的buffers项，c_ca_cyt,c_cam,c_trc,c_srm,c_slm用n时刻的值

    for i in range(0, NR):
        j_1 = -KPCAM * c_ca_cyt[i] * (BCAM - c_cam[i]) + KMCAM * c_cam[i]
        j_2 = -KPTRC * c_ca_cyt[i] * (BTRC - c_trc[i]) + KMTRC * c_trc[i]
        j_3 = -KPSRM * c_ca_cyt[i] * (BSRM - c_srm[i]) + KMSRM * c_srm[i]
        j_4 = -KPSLM * c_ca_cyt[i] * (BSLM - c_slm[i]) + KMSLM * c_slm[i]
        j_buffers[i] = j_1 + j_2 + j_3 + j_4
    return j_buffers


def cytosolic_buffers_equation(c_ca_cyt, c_cam, c_trc, c_srm, c_slm):
    # 更新[CaBi]

    for i in range(0, NR):
        # CCACYT用的上一步的值
        c_cam[i] = c_cam[i] + (KPCAM * c_ca_cyt[i] * (BCAM - c_cam[i]) - KMCAM * c_cam[i]) * DT
        c_trc[i] = c_trc[i] + (KPTRC * c_ca_cyt[i] * (BTRC - c_trc[i]) - KMTRC * c_trc[i]) * DT
        c_srm[i] = c_srm[i] + (KPSRM * c_ca_cyt[i] * (BSRM - c_srm[i]) - KMSRM * c_srm[i]) * DT
        c_slm[i] = c_slm[i] + (KPSLM * c_ca_cyt[i] * (BSLM - c_slm[i]) - KMSLM * c_slm[i]) * DT
    return c_cam, c_trc, c_srm, c_slm


def cytosolic_ca_equation(c_ca_cyt, j_dye, j_buffers, c_ca_store, KRYR2, iteration):
    # 更新[Ca2+]
    # 每个点计算完先放在med中
    # 这次迭代等所有点全部计算完再放到new中

    new_c_ca_cyt = np.zeros(NR, float)
    med_c_ca_cyt = np.zeros(NR, float)

    for i in range(0, NR):
        # 第一次迭代用n时刻的值
        new_c_ca_cyt[i] = c_ca_cyt[i]

    # 之后每次迭代用上一次迭代求的值
    for k in range(0, iteration):
        for i in range(0, NR):
            coeff = 3 * DCACYT * DT / DR

            # 第一个点
            if i == 0:
                r3 = (RADIUS[0] + RADIUS[1]) / 2.0
                r1 = 0
                fr3 = new_c_ca_cyt[1]

                # up:分子 down:分母
                j_diff_up = r3 * r3 * fr3 / (r3 ** 3 - r1 ** 3)
                j_diff_down = r3 * r3 / (r3 ** 3 - r1 ** 3)
                j_ryr_up = 3 * DT * KRYR2 * c_ca_store / (4 * pi * r3 ** 3)
                j_ryr_down = 3 * DT * KRYR2 / (4 * pi * r3 ** 3)

            # 最后一个点
            elif i == NR - 1:
                r3 = RADIUS[NR - 1]
                r1 = (RADIUS[NR - 1] + RADIUS[NR - 2]) / 2.0
                fr1 = new_c_ca_cyt[NR - 2]

                j_diff_up = r1 * r1 * fr1 / (r3 ** 3 - r1 ** 3)
                j_diff_down = r1 * r1 / (r3 ** 3 - r1 ** 3)
                j_ryr_up = 0.0
                j_ryr_down = 0.0

            # 其他的点
            else:
                r3 = (RADIUS[i] + RADIUS[i + 1]) / 2.0
                r1 = (RADIUS[i - 1] + RADIUS[i]) / 2.0
                fr3 = new_c_ca_cyt[i + 1]
                fr1 = new_c_ca_cyt[i - 1]

                j_diff_up = (r3 * r3 * fr3 + r1 * r1 * fr1) / (r3 ** 3 - r1 ** 3)
                j_diff_down = (r3 * r3 + r1 * r1) / (r3 ** 3 - r1 ** 3)
                j_ryr_up = 0.0
                j_ryr_down = 0.0

            # 每次迭代求得的值
            med_c_ca_cyt[i] = (coeff * j_diff_up + DT * j_dye[i] + DT * j_buffers[i] + j_ryr_up + c_ca_cyt[i]) / (
                    1 + coeff * j_diff_down + j_ryr_down)
        #  全部点都结束后赋值
        for i in range(0, NR):
            new_c_ca_cyt[i] = med_c_ca_cyt[i]
    # 最后一次赋给CCACYT
    for i in range(0, NR):
        c_ca_cyt[i] = new_c_ca_cyt[i]

    return c_ca_cyt


def cytosolic_gn_equation(c_caf, j_dye, iteration):
    # 更新[CaF]
    new_c_caf = np.zeros(NR, float)
    med_c_caf = np.zeros(NR, float)

    # 第一次迭代用n时刻的值
    for i in range(0, NR):
        new_c_caf[i] = c_caf[i]

    for k in range(0, iteration):
        for i in range(0, NR):
            coeff = 3 * DCAF * DT / DR
            if i == 0:
                r3 = (RADIUS[0] + RADIUS[1]) / 2.0
                r1 = 0
                gr3 = new_c_caf[1]

                j_diff_up = (r3 * r3 * gr3) / (r3 ** 3 - r1 ** 3)
                j_diff_down = r3 * r3 / (r3 ** 3 - r1 ** 3)
            elif i == NR - 1:
                r3 = RADIUS[NR - 1]
                r1 = (RADIUS[NR - 1] + RADIUS[NR - 2]) / 2.0
                gr1 = new_c_caf[NR - 2]

                j_diff_up = r1 * r1 * gr1 / (r3 ** 3 - r1 ** 3)
                j_diff_down = r1 * r1 / (r3 ** 3 - r1 ** 3)
            else:
                r3 = (RADIUS[i] + RADIUS[i + 1]) / 2.0
                r1 = (RADIUS[i - 1] + RADIUS[i]) / 2.0
                gr3 = new_c_caf[i + 1]
                gr1 = new_c_caf[i - 1]

                j_diff_up = (r3 * r3 * gr3 + r1 * r1 * gr1) / (r3 ** 3 - r1 ** 3)
                j_diff_down = (r3 * r3 + r1 * r1) / (r3 ** 3 - r1 ** 3)

            med_c_caf[i] = (coeff * j_diff_up - DT * j_dye[i] + c_caf[i]) / (1 + coeff * j_diff_down)

        for i in range(0, NR):
            new_c_caf[i] = med_c_caf[i]
    for i in range(0, NR):
        c_caf[i] = new_c_caf[i]
    return c_caf
