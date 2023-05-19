import random
from math import log, sqrt

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib import font_manager

font_manager.fontManager.addfont('../../fonts/time-simsun.ttf')

plt.rcParams['font.sans-serif'] = 'Times New Roman + Simsun'
plt.rcParams['axes.unicode_minus'] = False

# 读取数据
CCAF = 1 / 14
csv_data_ryr4 = pd.read_csv("avg_gn_jsr_4.csv", header=None)
csv_data_ryr6 = pd.read_csv("avg_gn_jsr_6.csv", header=None)
csv_data_ryr8 = pd.read_csv("avg_gn_jsr_8.csv", header=None)

# 刚开始先定义一维
temp_random = csv_data_ryr6

# 定义权重
weight_data = {'RyR4': 50, 'RyR6': 150, 'RyR8': 50}


def random_weight(weight_data):
    total = sum(weight_data.values())  # 权重求和
    ra = random.uniform(0, total)  # 在0与权重和之前获取一个随机数
    curr_sum = 0
    ret = None
    # keys = weight_data.iterkeys()  # 使用Python2.x中的iterkeys
    keys = weight_data.keys()  # 使用Python3.x中的keys
    for k in keys:
        curr_sum += weight_data[k]  # 在遍历中，累加当前权重值
        if ra <= curr_sum:  # 当随机数<=当前权重和时，返回权重key
            ret = k
            break
    return ret


for i in range(100):
    # 0-40ms时间随机,所以对应到文件就是随机数*200
    random_arr = int((random.random()) * 200)
    csv_data = csv_data_ryr6
    random_fill = csv_data.shift(random_arr).fillna(CCAF)
    # 迭代拼接
    temp_random = pd.concat([temp_random, random_fill], axis=1)

# 开始画图
plt.figure()
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.grid(linestyle="--")
random_x = temp_random.index.values
for j in range(temp_random.shape[1]):
    random_y = temp_random.iloc[:, j]
    plt.plot(random_x, random_y)
# plt.xlim((0, 500))
# plt.xticks([0, 100, 200, 300, 400, 500], [0, 20, 40, 60, 80, 100], fontsize=13)
plt.xlabel("time(ms)", fontsize=14)
plt.ylabel("concentration", fontsize=14)
plt.yticks(fontsize=13)
plt.savefig("随机100次.jpg")
plt.show()

# 画平均的
# 取个平均数
random_mean = temp_random.mean(axis=1).values

plt.figure()
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.grid(linestyle="--")

x = temp_random.index.values
y = random_mean

# from scipy.interpolate import make_interp_spline
#
# x_smooth = np.linspace(x.min(), x.max(), 300)
# y_smooth = make_interp_spline(x, y)(x_smooth)
plt.plot(x, y)
plt.plot(csv_data_ryr6.index.values, csv_data_ryr6[0].values)


# plt.plot(x_smooth, y_smooth)

# plt.xlim((0, 500))
# plt.xticks([0, 100, 200, 300, 400, 500], [0, 20, 40, 60, 80, 100], fontsize=13)
plt.xlabel("steps(*100)", fontsize=13)
plt.ylabel("concentration", fontsize=14)
plt.yticks(fontsize=13)
plt.savefig("result.jpg")
plt.show()
