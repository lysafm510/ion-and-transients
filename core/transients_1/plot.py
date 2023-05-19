import random
from math import log, sqrt

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from test3 import random_test
from matplotlib import font_manager

font_manager.fontManager.addfont('../../fonts/time-simsun.ttf')

plt.rcParams['font.sans-serif'] = 'Times New Roman + Simsun'
plt.rcParams['axes.unicode_minus'] = False

CCAF = 1 / 14
csv_data = pd.read_csv("avg_gn_jsr.csv", header=None)

temp_gauss = csv_data
temp_random = csv_data
temp_custom = csv_data

# ***********  高斯  **********************************************
n = []
for i in range(100):
    h = random.gauss(100, 67)
    n.append(h)

c = np.array(n)
c = np.sort(c)
c = c[(c >= 0.0) & (c <= 200.0)]
e = c.astype(int)
print(e)

# for i in range(0, 100):
for i in range(len(e)):
    # a = int((random.random()) * 100)
    a = e[i]
    b = csv_data.shift(a).fillna(CCAF)
    temp_gauss = pd.concat([temp_gauss, b], axis=1)

plt.figure()
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.grid(linestyle="--")

x = temp_gauss.index.values
for j in range(temp_gauss.shape[1]):
    y = temp_gauss.iloc[:, j]
    plt.plot(x, y)
plt.xlim((0, 500))
plt.xticks([0, 100, 200, 300, 400, 500], [0, 20, 40, 60, 80, 100], fontsize=13)
plt.xlabel("time(ms)", fontsize=14)
plt.ylabel("concentration", fontsize=14)
plt.yticks(fontsize=13)
plt.savefig("高斯100次.jpg")
plt.show()

# **************  随机  ***********************************************
for i in range(100):
    random_arr = int((random.random()) * 200)
    random_fill = csv_data.shift(random_arr).fillna(CCAF)
    temp_random = pd.concat([temp_random, random_fill], axis=1)

plt.figure()
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.grid(linestyle="--")
random_x = temp_random.index.values
for j in range(temp_random.shape[1]):
    random_y = temp_random.iloc[:, j]
    plt.plot(random_x, random_y)
plt.xlim((0, 500))
plt.xticks([0, 100, 200, 300, 400, 500], [0, 20, 40, 60, 80, 100], fontsize=13)
plt.xlabel("time(ms)", fontsize=14)
plt.ylabel("concentration", fontsize=14)
plt.yticks(fontsize=13)
plt.savefig("随机100次.jpg")
plt.show()
# ********************************************

# ********************************************************************
gauss_mean = temp_gauss.mean(axis=1).values
random_mean = temp_random.mean(axis=1).values
# print(row_mean)
plt.figure()

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.grid(linestyle="--")

x = temp_gauss.index.values
plt.plot(x, csv_data.values)
plt.plot(x, random_mean)
plt.plot(x, gauss_mean)
plt.legend(labels=["原曲线", "无规则随机", "高斯分布"])
plt.xlim((0, 500))
plt.xticks([0, 100, 200, 300, 400, 500], [0, 20, 40, 60, 80, 100], fontsize=13)
plt.xlabel("time(ms)", fontsize=14)
plt.ylabel("concentration", fontsize=14)
plt.yticks(fontsize=13)
plt.savefig("100次_对比.jpg")
plt.show()
