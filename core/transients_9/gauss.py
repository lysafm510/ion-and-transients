import csv
import random
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import font_manager

font_manager.fontManager.addfont('../../fonts/time-simsun.ttf')

plt.rcParams['font.sans-serif'] = 'Times New Roman + Simsun'
plt.rcParams['axes.unicode_minus'] = False

# 读取数据
CCAF = 1 / 14
CCA = 1
csv_data_ryr4 = pd.read_csv("avg_gn_jsr_4.csv", header=None)
csv_data_ryr6 = pd.read_csv("avg_gn_jsr_6.csv", header=None)
csv_data_ryr8 = pd.read_csv("avg_gn_jsr_8.csv", header=None)

csv_data_fn_4 = pd.read_csv("avg_fn_jsr_4.csv", header=None)
csv_data_fn_6 = pd.read_csv("avg_fn_jsr_6.csv", header=None)
csv_data_fn_8 = pd.read_csv("avg_fn_jsr_8.csv", header=None)

i_ca_fsr_4 = pd.read_csv("i_ca_fsr_4.csv", header=None)
i_ca_fsr_6 = pd.read_csv("i_ca_fsr_6.csv", header=None)
i_ca_fsr_8 = pd.read_csv("i_ca_fsr_8.csv", header=None)

i_ca_rsr_4 = pd.read_csv("i_ca_rsr_4.csv", header=None)
i_ca_rsr_6 = pd.read_csv("i_ca_rsr_6.csv", header=None)
i_ca_rsr_8 = pd.read_csv("i_ca_rsr_8.csv", header=None)

def insert(csv_data):
    xx = csv_data.index.values
    yy = csv_data[0].values
    return xx, yy


# 按照权重随机取数据
def prepare_data():
    # 定义权重
    weight = {'RyR4': 50, 'RyR6': 150, 'RyR8': 90}
    # 权重求和
    total_weight = sum(weight.values())
    # 在1与权重和之间获取一个随机数
    random_weight = round(random.uniform(1, total_weight))
    curr_sum = 0
    res = None
    # keys = weight.iterkeys()  # 使用Python2.x中的iterkeys
    keys = weight.keys()  # 使用Python3.x中的keys
    for key in keys:
        curr_sum += weight[key]  # 在遍历中，累加当前权重值
        if random_weight <= curr_sum:  # 当随机数<=当前权重和时，返回权重key
            res = key
            break
    return res


# 方法一
def gauss():
    # 拼接的第一列
    time_matrix_gn = csv_data_ryr6
    time_matrix_fn = csv_data_fn_6
    time_matrix_icaf = i_ca_fsr_6
    time_matrix_icar = i_ca_rsr_6

    # 取随机延长的时间
    n = []
    for i in range(100):
        # random.gauss(平均值,标准差)
        n.append(random.gauss(100, 100))
    extend_time = np.sort(np.array(n))
    # 筛选数据
    extend_time_array = (extend_time[(extend_time >= 0.0) & (extend_time <= 200.0)]).astype(int)
    print(extend_time_array)
    for i in range(len(extend_time_array)):
        # a = int((random.random()) * 100)
        time = extend_time_array[i]
        res = prepare_data()
        csv_data_gn = 0
        csv_data_fn = 0
        csv_data_icaf = 0
        csv_data_icar = 0

        if res == 'RyR4':
            csv_data_gn = csv_data_ryr4
            csv_data_fn = csv_data_fn_4
            csv_data_icaf = i_ca_fsr_4
            csv_data_icar = i_ca_rsr_4
        elif res == 'RyR6':
            csv_data_gn = csv_data_ryr6
            csv_data_fn = csv_data_fn_6
            csv_data_icaf = i_ca_fsr_6
            csv_data_icar = i_ca_rsr_6
        elif res == 'RyR8':
            csv_data_gn = csv_data_ryr8
            csv_data_fn = csv_data_fn_8
            csv_data_icaf = i_ca_fsr_8
            csv_data_icar = i_ca_rsr_8

        time_matrix_gn = pd.concat([time_matrix_gn, csv_data_gn.shift(time).fillna(CCAF)], axis=1)
        time_matrix_fn = pd.concat([time_matrix_fn, csv_data_fn.shift(time).fillna(CCA)], axis=1)
        time_matrix_icaf = pd.concat([time_matrix_icaf, csv_data_icaf.shift(time).fillna(csv_data_icaf[0])], axis=1)
        time_matrix_icar = pd.concat([time_matrix_icar, csv_data_icar.shift(time).fillna(csv_data_icar[0])], axis=1)

    return time_matrix_gn, time_matrix_fn, time_matrix_icaf, time_matrix_icar


# 方法二
def custom():
    time_matrix = csv_data_ryr6
    a = np.linspace(0, 80, 100)
    b = np.linspace(80, 120, 100)
    c = np.linspace(120, 200, 10)
    d = np.concatenate((a[:-1], b[:-1], c[:-1]))
    # print(d)
    for i in range(len(d)):
        # 随机数偏移量,并填充
        res = prepare_data()
        csv_data = 0
        if res == 'RyR4':
            csv_data = csv_data_ryr4
        elif res == 'RyR6':
            csv_data = csv_data_ryr6
        elif res == 'RyR8':
            csv_data = csv_data_ryr8
        # 迭代拼接
        time_matrix = pd.concat([time_matrix, csv_data.shift(int(d[i])).fillna(CCAF)], axis=1)
    return time_matrix


time_matrix_gn, time_matrix_fn, time_matrix_icaf, time_matrix_icar = gauss()

gn_x = time_matrix_gn.index.values
# 求均值
gn_y = time_matrix_gn.mean(axis=1).values

fn_x = time_matrix_fn.index.values
# 求均值
fn_y = time_matrix_fn.mean(axis=1).values

icaf_x = time_matrix_icaf.index.values
# 求均值
icaf_y = time_matrix_icaf.mean(axis=1).values

icar_x = time_matrix_icar.index.values
# 求均值
icar_y = time_matrix_icar.mean(axis=1).values

np.savetxt('avg_gn_jsr_gauss_new_2.csv', gn_y, fmt='%.15f')
np.savetxt('avg_fn_jsr_gauss_new_2.csv', fn_y, fmt='%.15f')
np.savetxt('i_ca_fsr_gauss_new_2.csv', icaf_y, fmt='%.15f')
np.savetxt('i_ca_rsr_gauss_new_2.csv', icar_y, fmt='%.15f')

print(type(gn_y))
plt.figure()
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.grid(linestyle="--")
plt.plot(csv_data_ryr4.index.values, csv_data_ryr4[0].values, linewidth=1)
plt.plot(csv_data_ryr6.index.values, csv_data_ryr6[0].values, linewidth=1)
plt.plot(csv_data_ryr8.index.values, csv_data_ryr8[0].values, linewidth=1)
plt.plot(gn_x, gn_y, linewidth=1)
plt.title("[CaF]")
plt.ylabel("concentration", fontsize=14)
plt.xlabel("time(ms)", fontsize=14)
plt.xlim((0, 900))
plt.xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900],
           [0, 20, 40, 60, 80, 100, 120, 140, 160, 180], fontsize=13)
plt.yticks(fontsize=13)
plt.savefig("gauss_new.jpg")
plt.show()
