import random
from math import log, sqrt

import numpy as np
from matplotlib import pyplot as plt


# x = []
# i = 0
# while i < 20:
#     i = i + 1
#     x1 = random.randint(-20, 20)  # 随机生成整数,其范围在(-20,20)区间
#     x.append(x1)
# # 设立y轴的值相等就可以绘制一维图了
# y = []
# for i in range(20):
#     y.append(0)
# plt.scatter(x, y, edgecolors='red')
# plt.xlim(-20, 20)
# plt.ylim(-1, 1)
# plt.show()


def random_test():
    u = sqrt(-2 * log(1 / 3))
    print(u)

    j = np.random.randn(100)
    c = j[(j <= u) & (j >= -u)]
    print(len(c))

    d = c * 100 / u + 100
    print(d)
    e = d.astype(int)
    print(e)

    # plot(e)
    return e


# def plot(e):
#     y = []
#     for i in range(len(e)):
#         y.append(e[i])
#     plt.scatter(e, y)
#     plt.xlim(-1000, 1000)
#     plt.ylim(-1000,1000)
#     plt.show()
