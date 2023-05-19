from math import sqrt
import numpy as np
from blink_coeff_matrix import LN_MATRIX, CTRL_AREA, NMAX, TRIANGLE_NUMBER, ABC_MATRIX, AREA, TOTAL_AREA
from constant import B_INNER, B_INFLOW, B_OUTFLOW, K2, K1, DCAJSR, F, BCSQ, KDCSQ, DCAFSR, CCAFSR, DCAF, \
    DT, H_JSR, UNITEC, MOLNUM
from grid_info import NOD, NPOCH, NP, PX, PY, NE


def point_type(num_of_tri):
    """
    判断三角形不同点的类型并返回数组
    :param num_of_tri: 三角形序号
    :return: inner_point, out_boundary, in_boundary, in_L, out_L
    """
    inner_point = []
    out_boundary = []
    in_boundary = []
    in_l = 0.
    out_l = 0.
    for i in range(0, 3):
        n = NOD[i, num_of_tri]
        if NPOCH[n] == B_INNER:
            inner_point.append(n)  # 存的点位置
        elif NPOCH[n] == B_INFLOW:
            in_boundary.append(n)
        elif NPOCH[n] == B_OUTFLOW:
            out_boundary.append(n)

    if len(in_boundary) == 2:
        for i in range(0, 3):
            n = NOD[i, num_of_tri]
            if (n != in_boundary[0]) and (n != in_boundary[1]):
                in_l = LN_MATRIX[num_of_tri][i][0]

    if len(out_boundary) == 2:
        for i in range(0, 3):
            n = NOD[i, num_of_tri]
            if (n != out_boundary[0]) & (n != out_boundary[1]):
                out_l = LN_MATRIX[num_of_tri][i][0]

    return inner_point, out_boundary, in_boundary, in_l, out_l


def sort(n1, n2, n3, i):
    """
    :param n1:
    :param n2:
    :param n3:
    :param i: 中心点
    :return:Nod_2, Nod_3, Column_2, Column_3
    """
    nod_2 = -1
    nod_3 = -1
    column_2 = -1
    column_3 = -1
    if i == n1:
        nod_2, nod_3 = n2, n3
        column_2, column_3 = 1, 2
    elif i == n2:
        nod_2, nod_3 = n3, n1
        column_2, column_3 = 2, 0
    elif i == n3:
        nod_2, nod_3 = n1, n2
        column_2, column_3 = 0, 1
    return nod_2, nod_3, column_2, column_3


def blink_fn_equation(fn, gn, CCACYTO, DCARYR, iteration):
    """
    构建钙离子 n+1 步的方程
    """

    new_fn = np.zeros(NP, float)  # 存每一次迭代的结果fn
    med_fn = np.zeros(NP, float)  # 临时变量
    for k in range(0, iteration):
        for i in range(0, NP):  # i代表点， j第j个三角形 N1代表该三角形的其中一点
            cal_one = 0.
            cal_two = (K2 * gn[i] - K1 * fn[i] * (F - gn[i])) * CTRL_AREA[i]  # 系数
            cal_three = 0.
            cal_four = 0.

            cal_five = 0.  # 分子入流出流
            cal_six = 0.  # 分母入流出流

            for j in range(0, NMAX):
                if TRIANGLE_NUMBER[i][j][0] != -1:
                    tri_pos = TRIANGLE_NUMBER[i][j][0]  # 点i的相关的 第j个三角形的绝对序号
                    nod_index = TRIANGLE_NUMBER[i][j][1]  # 点i在NOD数组中这个三角形中相对第几个点

                    # i nod_index     j tri_pos
                    # tri_pos 为j三角形绝对序号
                    # nod_index 为中心点绝对位置i在每一个三角形中对应的相对序号，0或1或2
                    n1 = NOD[0, tri_pos]  # 点i的相关的第j个三角形的每个点
                    n2 = NOD[1, tri_pos]
                    n3 = NOD[2, tri_pos]

                    nod_2, nod_3, column_2, column_3 = sort(n1, n2, n3, i)

                    if k == 0:
                        cca1_nod1 = fn[i]
                        cca1_nod2 = fn[nod_2]
                        cca1_nod3 = fn[nod_3]
                    else:
                        cca1_nod1 = (fn[i] + new_fn[i]) / 2.0
                        cca1_nod2 = (fn[nod_2] + new_fn[nod_2]) / 2.0
                        cca1_nod3 = (fn[nod_3] + new_fn[nod_3]) / 2.0
                    # 求的不带边界点的
                    # ABC_MATRIX : 第tri_pos个三角形，相对第(0，1，2)个点，第3维度(1 a;2 b;3 c)
                    # LN_MATRIX : 第tri_pos个三角形，相对第(0，1，2)个点，第3维度(1 L;2 Nix;3 Niy)
                    # nod_index 中心点   column_2、column_3 另外两个点    nod_2、nod_3  另外两个点绝对位置
                    # cal_one 从0开始

                    # nod_2,nod_3只要有一个是内部结点，就计算这个三角形
                    # nod_2,nod_3都是边界点（>0），跳过这个三角形
                    if NPOCH[nod_2] == B_INNER or NPOCH[nod_3] == B_INNER:
                        cal_one = cal_one + DCAJSR * (LN_MATRIX[tri_pos][nod_index][1] *
                                                      (cca1_nod2 * ABC_MATRIX[tri_pos][column_2][1] +
                                                       cca1_nod3 * ABC_MATRIX[tri_pos][column_3][1])
                                                      + LN_MATRIX[tri_pos][nod_index][2] *
                                                      (cca1_nod2 * ABC_MATRIX[tri_pos][column_2][2] +
                                                       cca1_nod3 * ABC_MATRIX[tri_pos][column_3][2])
                                                      ) * LN_MATRIX[tri_pos][nod_index][0]

                    # 计算分子分母公共项
                    avg_cca = (cca1_nod1 + cca1_nod2 + cca1_nod3) / 3.0
                    cal_three = cal_three + (1 + BCSQ * KDCSQ / ((KDCSQ + avg_cca) ** 2)) * AREA[tri_pos]

                    # 判断点的类型,输入第tri_pos个三角形，第i个点
                    inner_point, out_boundary_point, in_boundary_point, in_l, out_l = point_type(tri_pos)

                    if k == 0:
                        cca2_nod2 = fn[nod_2]
                        cca2_nod3 = fn[nod_3]
                    else:
                        cca2_nod2 = new_fn[nod_2]
                        cca2_nod3 = new_fn[nod_3]

                    if (len(out_boundary_point)) == 2:  # 为出流边界点时
                        fk = (cca2_nod2 + cca2_nod3) / 3.0
                        cal_five = cal_five + DCARYR * (CCACYTO - fk) * out_l  # 分子
                        cal_six = cal_six + (DCARYR * out_l) / 3.0  # 分母

                    if (len(in_boundary_point)) == 2:  # 为入流边界点时
                        fj = (cca2_nod2 + cca2_nod3) / 3.0
                        cal_five = cal_five + DCAFSR * (CCAFSR - fj) * in_l  # 分子
                        cal_six = cal_six + (DCAFSR * in_l) / 3.0  # 分母

                    # 计算分子分母公共项
                    # nod_2,nod_3只要有一个是内部结点，就计算这个三角形
                    # nod_2,nod_3都是边界点（>0），跳过这个三角形
                    if NPOCH[nod_2] == B_INNER or NPOCH[nod_3] == B_INNER:
                        cal_four = cal_four + DCAJSR * (1 / 2) * (
                                LN_MATRIX[tri_pos][nod_index][1] * ABC_MATRIX[tri_pos][nod_index][1] +
                                LN_MATRIX[tri_pos][nod_index][2] * ABC_MATRIX[tri_pos][nod_index][2]) * \
                                   LN_MATRIX[tri_pos][nod_index][0]
                else:
                    break
            # 中间结果
            med_fn[i] = (DT * cal_one + DT * cal_two + fn[i] * cal_three + DT * fn[i]
                         * cal_four + DT * cal_five) / (cal_three - DT * cal_four + DT * cal_six)
        # 当每次迭代所有点都计算完后，再赋值
        for i in range(0, NP):
            new_fn[i] = med_fn[i]
    return new_fn


def blink_gn_equation(fn, gn, iteration):
    """
    构建Gn方程
    """
    new_gn = np.zeros(NP, float)
    med_gn = np.zeros(NP, float)

    for k in range(0, iteration):
        for i in range(0, NP):  # i代表点，m代表该点三角形的个数 j第j个三角形 N1代表该三角形的其中一点
            cal_one = 0.
            cal_two = (K1 * fn[i] * (F - gn[i]) - K2 * gn[i]) * CTRL_AREA[i]
            cal_three = 0.

            for j in range(0, NMAX):
                if TRIANGLE_NUMBER[i][j][0] != -1:
                    tri_pos = TRIANGLE_NUMBER[i][j][0]
                    nod_index = TRIANGLE_NUMBER[i][j][1]

                    n1 = NOD[0, tri_pos]
                    n2 = NOD[1, tri_pos]
                    n3 = NOD[2, tri_pos]

                    nod_2, nod_3, column_2, column_3 = sort(n1, n2, n3, i)

                    if k == 0:
                        gn_nod2 = gn[nod_2]
                        gn_nod3 = gn[nod_3]
                    else:
                        gn_nod2 = (gn[nod_2] + new_gn[nod_2]) / 2.0
                        gn_nod3 = (gn[nod_3] + new_gn[nod_3]) / 2.0
                    # 求的不带边界点的
                    # abcMatrix : 第tri_pos个三角形，相对第(0，1，2)个点，第3维度(1 a;2 b;3 c)
                    # nlMatrix : 第tri_pos个三角形，相对第(0，1，2)个点，第3维度(1 L;2 Nix;3 Niy)
                    # nod_index 中心点   column_2、column_3 另外两个点    nod_2、nod_3  另外两个点绝对位置

                    # nod_2,nod_3只要有一个是内部结点，就计算这个三角形
                    # nod_2,nod_3都是边界点（>0），跳过这个三角形
                    if NPOCH[nod_2] == B_INNER or NPOCH[nod_3] == B_INNER:
                        cal_one = cal_one + DCAF * (LN_MATRIX[tri_pos][nod_index][1] *
                                                    (gn_nod2 * ABC_MATRIX[tri_pos][column_2][1] +
                                                     gn_nod3 * ABC_MATRIX[tri_pos][column_3][1])
                                                    + LN_MATRIX[tri_pos][nod_index][2] *
                                                    (gn_nod2 * ABC_MATRIX[tri_pos][column_2][2] +
                                                     gn_nod3 * ABC_MATRIX[tri_pos][column_3][2])
                                                    ) * LN_MATRIX[tri_pos][nod_index][0]
                        cal_three = cal_three + 0.5 * DCAF * (
                                LN_MATRIX[tri_pos][nod_index][1] * ABC_MATRIX[tri_pos][nod_index][1] +
                                LN_MATRIX[tri_pos][nod_index][2] * ABC_MATRIX[tri_pos][nod_index][2]) * \
                                    LN_MATRIX[tri_pos][nod_index][0]

                else:
                    break
            med_gn[i] = (DT * cal_one + DT * cal_two + gn[i] * CTRL_AREA[i] + gn[i] * DT * cal_three) / (
                    CTRL_AREA[i] - DT * cal_three)
        for i in range(0, NP):
            new_gn[i] = med_gn[i]
    for i in range(0, NP):
        gn[i] = new_gn[i]

    return gn


def average_concentration(conc_matrix):
    """
    求平均钙离子/荧光钙浓度
    :param conc_matrix: 钙离子/荧光钙数值
    :return: avg_conc 平均浓度
    """
    total_conc = 0.0
    avg_conc = 0.0
    for i in range(0, NP):
        total_conc = total_conc + conc_matrix[i] * CTRL_AREA[i]
    total_conc = total_conc / 3
    if TOTAL_AREA > 0.00000000000001:
        avg_conc = total_conc / TOTAL_AREA
    return avg_conc


def ca_current(fn, DCARYR):
    """
    求钙离子电流
    """
    ca_in = 0.0
    ca_out = 0.0
    for i in range(0, NE):
        for j in range(0, 3):
            p1 = NOD[(j + 1) % 3, i]
            p2 = NOD[(j + 2) % 3, i]
            #  out current
            if (NPOCH[p1] == B_OUTFLOW) and (NPOCH[p2] == B_OUTFLOW):
                length = sqrt((PX[p1] - PX[p2]) ** 2 + (PY[p1] - PY[p2]) ** 2)
                ca_out = ca_out + length * (fn[p1] + fn[p2]) / 2
            #  in current
            if NPOCH[p1] == B_INFLOW and NPOCH[p2] == B_INFLOW:
                length = sqrt((PX[p1] - PX[p2]) ** 2 + (PY[p1] - PY[p2]) ** 2)
                ca_in = ca_in + length * (CCAFSR - (fn[p1] + fn[p2]) / 2)
    ca_out = ca_out * DCARYR * H_JSR
    ca_in = ca_in * DCAFSR * H_JSR
    i_ca_ryr = ca_out * UNITEC * 2 * MOLNUM * (10.0 ** -11)
    i_ca_fsr = ca_in * UNITEC * 2 * MOLNUM * (10.0 ** -11)

    return i_ca_ryr, i_ca_fsr


def cca_ryr_inside_store(fn):
    ca_out = 0.0
    arc_len = 0.0
    for i in range(0, NE):
        for j in range(0, 3):
            p1 = NOD[(j + 1) % 3, i]
            p2 = NOD[(j + 2) % 3, i]
            #  out current
            if (NPOCH[p1] == B_OUTFLOW) and (NPOCH[p2] == B_OUTFLOW):
                length = sqrt((PX[p1] - PX[p2]) ** 2 + (PY[p1] - PY[p2]) ** 2)
                ca_out = ca_out + length * (fn[p1] + fn[p2]) / 2
                arc_len = arc_len + length
    cca_store = ca_out / arc_len

    return cca_store
