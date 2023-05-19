import os
import numpy as np
from math import sqrt

from constant import TRIANGLE_PATH
from grid_info import NP, NE, NOD, PX, PY


def cal_area():
    """
    计算每个三角形面积/每个点控制面积
    """
    ctrl_area = np.zeros(NP, float)  # 每一个点的控制面积
    total_area = 0.0  # 总面积
    det = np.zeros(NE, float)  # 存储每个三角形的D
    area = np.zeros(NE, float)  # 存储每个三角形的A
    matrix = np.zeros((3, 3))  # 矩阵用来计算行列式
    count = np.zeros(NP, float)  # 计算每个点相关的三角形的个数

    for e in range(0, NE):  # NE为三角形单元的数量,计算D和V
        n1 = NOD[0, e]
        n2 = NOD[1, e]
        n3 = NOD[2, e]
        count[n1] = count[n1] + 1
        count[n2] = count[n2] + 1
        count[n3] = count[n3] + 1  # 计算每个点相关的三角形的个数
        for i in range(0, 3):
            matrix[i, 0] = 1.0
            matrix[i, 1] = PX[NOD[i, e]]
            matrix[i, 2] = PY[NOD[i, e]]
        det1 = np.linalg.det(matrix)  # 计算行列式
        a = det1 / 2.0  # 计算三角形面积
        det[e] = det1
        area[e] = a
        total_area = total_area + a

        for i in range(0, 3):
            ctrl_area[NOD[i, e]] = ctrl_area[NOD[i, e]] + a  # 每一点的控制面积
    nmax = int(max(count))
    return ctrl_area, total_area, det, area, nmax


def cal_abc_l_n(det):
    """
    计算a,b,c,N,L系数矩阵
    :param det: 存储每个三角形的D
    :return abcMatrix, nlMatrix
    """
    abc_matrix = np.zeros((NE, 3, 3))  # 存储每个三角形的a,b,c系数
    ln_matrix = np.zeros((NE, 3, 3))  # 存储每个三角形的l,nix,niy系数
    for i in range(0, NE):
        for j in range(0, 3):
            p2 = NOD[(j + 1) % 3, i]
            p3 = NOD[(j + 2) % 3, i]
            d = det[i]

            # 第i个三角形，每个点为中心的a,b,c
            abc_matrix[i, j, 0] = (PX[p2] * PY[p3] - PX[p3] * PY[p2]) / d
            abc_matrix[i, j, 1] = (PY[p2] - PY[p3]) / d
            abc_matrix[i, j, 2] = (PX[p3] - PX[p2]) / d

            # 第i个三角形，每个点为中心的L，Nix,Niy
            ln_matrix[i, j, 0] = sqrt((PX[p2] - PX[p3]) ** 2 + (PY[p2] - PY[p3]) ** 2)
            ln_matrix[i, j, 1] = (PY[p3] - PY[p2]) / ln_matrix[i, j, 0]
            ln_matrix[i, j, 2] = (PX[p2] - PX[p3]) / ln_matrix[i, j, 0]

    return abc_matrix, ln_matrix


def not_empty(s):
    return s and s.strip()


def search_triangle(nmax):
    """
    存储某一点相邻的所有三角形单元
    """
    triangle_number = np.full([NP, nmax, 2], -1)  # 存放围绕此点的三角形单元编号
    if not os.path.exists("../../parameters/triangle_number/" + TRIANGLE_PATH):
        print("write...triangle number...")
        os.makedirs("../../parameters/triangle_number/" + TRIANGLE_PATH)
    else:
        print("load...triangle number...")
    for n in range(0, NP):  # 遍历所有的点
        index = 0
        # 先检查有没有文件
        # 1.没有文件->写
        if not os.path.exists("../../parameters/triangle_number/" + TRIANGLE_PATH + "/POINT" + str(n + 1) + ".dat"):

            with open("../../parameters/triangle_number/" + TRIANGLE_PATH + "/POINT" + str(n + 1) + ".dat", "w") as file_object:
                for E in range(0, NE):  # 遍历每个三角形
                    for i in range(0, 3):  # 遍历三角形的三个顶点
                        if n == NOD[i][E]:
                            triangle_number[n][index][0] = E
                            triangle_number[n][index][1] = i
                            file_object.write(str(E) + " ")
                            file_object.write(str(i) + "\n")
                            index = index + 1
                            break
            file_object.close()
        # 2.有文件->直接读取
        else:
            with open("../../parameters/triangle_number/" + TRIANGLE_PATH + "/POINT" + str(n + 1) + ".dat", "r") as file_object:
                for line in file_object.readlines():
                    current_line = list(filter(not_empty, line.strip("\n").split(" ")))
                    triangle_number[n][index][0] = current_line[0]
                    triangle_number[n][index][1] = current_line[1]
                    index = index + 1
            file_object.close()
    print("write/load...triangle number...over")
    return triangle_number


CTRL_AREA, TOTAL_AREA, DET, AREA, NMAX = cal_area()
ABC_MATRIX, LN_MATRIX = cal_abc_l_n(DET)
TRIANGLE_NUMBER = search_triangle(NMAX)
