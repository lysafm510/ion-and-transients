import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from numpy.linalg import det

global grid, PX, PY, NP, nod, NE, NOD

# 读取文件夹
data_path = "test10_测试gn上升_带电荷力"
# fn或者gn
type = "gn"
# 读取网格路径
grid_path = "beibei"
# 中间参数保存路径
param_path = "beibei"


# ***************************************************************

def grid_info(grid_path):
    global grid, PX, PY, NP, nod, NE, NOD

    grid = np.loadtxt("../parameters/grid/blink/" + grid_path + "/gridt.dat", dtype=np.float64)
    PX = grid[:, 0]
    PY = grid[:, 1]
    NP = len(grid)

    nod = np.loadtxt("../parameters/grid/blink/" + grid_path + "/nod.dat", dtype=np.int64)
    NE = len(nod)
    NOD = nod.T - 1


def preparation():
    """
    准备工作
    """
    prepare = np.zeros([300, 4], float)  # 第一列存x坐标，第二列存y坐标，第三列存在哪个三角形中，第四列存该三角形面积

    arr1 = np.zeros([3, 3], float)  # 临时数组
    arr2 = np.zeros([3, 3], float)  # 临时数组
    area = np.zeros(3, float)  # 求三个三角形面积
    num = 0

    for i in range(0, 300):  # 300个点，-299开始
        prepare[i][0] = -299 + num * 2
        prepare[i][1] = 0
        # 依次判断 i 点在哪个网格三角形中
        for j in range(0, NE):
            # 计算该网格三角形的面积
            for k in range(0, 3):
                arr1[k, 0] = 1
                arr1[k, 1] = PX[NOD[k, j]]
                arr1[k, 2] = PY[NOD[k, j]]
            total_area = abs(det(arr1))

            # 分别计算点插入后，三个小三角形的面积，求和。三个小三角形面积和=大面积，说明在这个点在这个三角形中
            for k in range(0, 3):
                arr2[0, 0] = 1
                arr2[0, 1] = PX[NOD[(k + 1) % 3, j]]
                arr2[0, 2] = PY[NOD[(k + 1) % 3, j]]
                arr2[1, 0] = 1
                arr2[1, 1] = PX[NOD[(k + 2) % 3, j]]
                arr2[1, 2] = PY[NOD[(k + 2) % 3, j]]
                arr2[2, 0] = 1
                arr2[2, 1] = prepare[i][0]
                arr2[2, 2] = prepare[i][1]
                area[k] = abs(det(arr2))
            # 三个小三角形的面积
            sum_area = area[0] + area[1] + area[2]
            # 差值在一个很小的范围内
            if abs(total_area - sum_area) < 10e-8:
                prepare[i][2] = j
                prepare[i][3] = total_area
                break  # 找到后跳出循环
        num = num + 1
        print("[" + str(int(prepare[i][0])) + "," + str(int(prepare[i][1])) + "]点计算完毕")

    data = pd.DataFrame(prepare)
    data.to_csv("../parameters/plot/" + param_path + ".csv", index=False, header=False)

    return prepare


def travers(data_path, type, param_path):
    """
    遍历
    """
    if not os.path.exists("../parameters/plot/" + param_path + ".csv"):
        prepare = np.array(preparation())
    else:
        prepare = np.array(pd.read_csv("../parameters/plot/" + param_path + ".csv", header=None))
    arr3 = np.zeros([3, 3], float)  # 临时数组
    file_list = os.listdir("../data/" + data_path + "/blink/" + type)
    result = np.zeros((len(prepare), len(file_list)))
    # 遍历每个文件
    for i in range(0, len(file_list)):
        df = pd.read_csv("../data/" + data_path + "/blink/" + type + "/" + file_list[i], dtype=str, header=None)
        concentration = 0
        if type == "fn":
            concentration = np.array(df.iloc[:len(df) - 4, 0]).astype('float64')
        elif type == "gn":
            concentration = np.array(df.iloc[:len(df) - 2, 0]).astype('float64')
        # 遍历每个点
        for j in range(0, len(prepare)):
            x = int(prepare[j][0])
            y = int(prepare[j][1])
            tri = int(prepare[j][2])
            area = prepare[j][3]
            value = 0  # 存这个点的浓度
            for k in range(0, 3):
                temp = concentration[NOD[k, tri]]
                arr3[0, 0] = 1
                arr3[0, 1] = PX[NOD[(k + 1) % 3, tri]]
                arr3[0, 2] = PY[NOD[(k + 1) % 3, tri]]
                arr3[1, 0] = 1
                arr3[1, 1] = PX[NOD[(k + 2) % 3, tri]]
                arr3[1, 2] = PY[NOD[(k + 2) % 3, tri]]
                arr3[2, 0] = 1
                arr3[2, 1] = x
                arr3[2, 2] = y
                value = value + abs(det(arr3)) * temp / area

            result[j][i] = value
    return result


def plot(data_path, data):
    fig, ax = plt.subplots()
    # imshow 绘制热图
    im = ax.imshow(data)
    plt.colorbar(im)

    plt.yticks([300, 250, 200, 150, 100, 50, 0], [-300, -200, -100, 0, 100, 200, 300])
    plt.xticks([0, 100, 200, 300, 400, 500], [0, 0.02, 0.04, 0.06, 0.08, 0.1])
    plt.xlabel("time(s)")
    plt.ylabel("x")
    plt.savefig("../figure/square/" + data_path + ".jpg")
    plt.show()


if __name__ == '__main__':
    grid_info(grid_path)
    plot(data_path, travers(data_path, type, param_path))
