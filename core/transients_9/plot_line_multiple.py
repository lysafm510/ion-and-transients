import matplotlib.pyplot as plt
import csv
from matplotlib import font_manager
import time

font_manager.fontManager.addfont('../../fonts/time-simsun.ttf')

plt.rcParams['font.sans-serif'] = 'Times New Roman + Simsun'
# plt.rcParams['font.sans-serif'] = 'Times New Roman'
plt.rcParams['axes.unicode_minus'] = False
# plt.rcParams["font.family"] = "Microsoft Yahei"
# plt.rcParams['font.sans-serif'] = ['SimHei']

# path1 = "20220517_random4_DCAFSR_33"
# path2 = "20220514_random6_DCAFSR_33"
# path3 = "20220517_random8_DCAFSR_33"

path1 = "avg_gn_jsr_4.csv"
path2 = "avg_gn_jsr_6.csv"
path3 = "avg_gn_jsr_8.csv"
path4 = "avg_gn_jsr_gauss_new_2.csv"

# path1 = "20220607_random4_ISO"
# path2 = "20220607_random6_ISO"
# path3 = "20220607_random8_ISO"

# path1 = "20220401_4RyR_DCARYR_25%"
# path2 = "20220401_4RyR_DCARYR_50%"
# path3 = "20220401_4RyR_DCARYR_75%"
# path4 = "20220401_4RyR_DCARYR_100%"

# path1 = "20220514_random6_DCAFSR_20"
# path2 = "20220514_random6_DCAFSR_25"
# path3 = "20220514_random6_DCAFSR_33"
# path4 = "20220514_random6_DCAFSR_50"
# path5 = "20220514_random6_DCAFSR_100"

# path1 = "20220501_center4"
# path2 = "20220501_random4"

# path1 = "20220709_caffeine_0"
# path2 = "20220709_caffeine_33"
# path3 = "20220623_random25_DCAFSR_0"
# path4 = "20220623_random25_DCAFSR_33"

# path1 = "20220717_DCAJSR_0.01"
# path2 = "20220704_center4_DCAJSR_0.1"
# path3 = "20220709_DCAJSR_0.2"
# path4 = "20220709_DCAJSR_0.5"
# path5 = "20220709_DCAJSR_1"
# path4 = "20220709_DCAJSR_1"
# path5 = "20220709_DCAJSR_2"
# path6 = "20220709_DCAJSR_5"
# path7 = "20220704_center4_DCAJSR_10"

# path1 = "20220704_center4_K1K2_0.1"
# path2 = "20220709_center4_DCAF_0.1"
# path3 = "20220709_center4_K1K2_DCAF_0.1"
# path4 = "20220713_center4_DCAF_0.1"
# path5 = "20220713_center4_K1K2_DCAF_0.1"

# path1 = "20220717_DCAF_6_0.1"
# path2 = "20220717_DCAF_6_0.01"
# path3 = "20220717_K1K2_0.01"
# path4 = "20220717_K1K2_0.001"
# path5 = "20220717_K1K2_0.01_DCAF_0.01"
# path1 = "20220709_DCAJSR_1"
# path2 = "20220801_center4_DCAF_0.1"
# path3 = "20220801_center4_DCAF_0.01"
# path4 = "20220801_center4_DCAF_0.001"
# path5 = "20220801_center4_DCAF_0.0001"
# path6 = "20220801_center4_DCAF_0.00001"
# path7 = "20220801_center4_DCAF_0.000001"

# path1 = "20220709_DCAJSR_1"
# path2 = "20220704_center4_K1K2_0.1"
# path3 = "20220801_center4_K1K2_0.05"
# path4 = "20220801_center4_K1K2_0.02"
# path5 = "20220801_center4_K1K2_0.01"
# path6 = "20220811_center4_K1K2_0.005"
# path7 = "20220814_center4_K1K2_0.003"
# path8 = "20220814_center4_K1K2_0.002"
# path9 = "20220801_center4_K1K2_0.001"

# path1 = "20220704_center2_DCARYR_2"
# path2 = "20220501_center4"

# path1 = "20220709_caffeine_0"
# path2 = "20220709_caffeine_33"

# path1 = "20220717_DCAF_6_0.01"
# path2 = "20220717_K1K2_0.01"
# path3 = "20220717_K1K2_0.01_DCAF_0.01"

# paths = [path1, path2, path3, path4, path5, path6, path7, path8, path9]
paths = [path1, path2, path3, path4]
# paths = [path1, path2]

# label_list = ["0.1", "0.2", "0.5", "1", "2", "5", "10", "old"]
# label_list = ["2*10^7", "2*10^6", "2*10^5", "2*10^4", "2*10^3", "2*10^2", "2*10^1"]
# label_list = ["1", "0.1", "0.05", "0.02", "0.01", "0.005", "0.003", "0.002", "0.001"]
# label_list = ["DCAF = 2*10^4", "K1 = 488, K2 = 195.2", "(1)&(2)","a"]
# label_list = ["2RyRs , DCARYR*0.5", "4RyRs"]

# label_list = ["25 %", "50 %", "75 %", "100 %"]
# label_list = ["center 4 RyRs", "random 4 RyRs"]
# label_list = ["20 %", "25 %", "33 %", "50 %", "100 %"]
label_list = ["4 RyRs", "6 RyRs", "8 RyRs", "curve fitting"]

# 文件描述
filename = ""
title = ""


# *****************************************************************************
def plot(label_list, type, filename, title, paths):
    """
    绘制x:time y:concentration曲线图

    :param label_list: legend图例
    :param type: ['fn','gn']
    :param filename: 文件名前缀
    :param title: 标题前缀
    :param paths: 读取数据路径数组
    :return:
    """
    n = len(paths)
    plt.figure(figsize=plt.figaspect(3 / 4))
    max = 0
    save_steps = 100
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.grid(linestyle="--")
    for i in range(0, n):
        x = []
        y = []
        j = 0

        with open(paths[i], 'r') as ca_csvfile1:
            plots = csv.reader(ca_csvfile1, delimiter='\t')
            for row in plots:
                x.append(j * save_steps)
                y.append(float(row[0]))
                j = j + 1

                # if j * save_steps > 10000:
                #     break

        plt.plot(x, y, linewidth=1)

        if j * save_steps > max:
            max = j * save_steps

    # parameter = "DCAFSR"
    if type == "fn":
        # figure_title = "[Ca] - adjust parameter \"" + parameter + "\""
        # figure_title = "[Ca] - adjust 4RyRs positions"
        figure_title = "[Ca] - adjust the number of RyRs"
    else:
        # figure_title = "[CaF] - adjust parameter \"" + parameter + "\""
        # figure_title = "[CaF] - adjust 4RyRs positions"
        figure_title = "[CaF] - curve fitting"
    plt.title(figure_title, fontsize=14)
    plt.legend(labels=label_list, fontsize=12, loc=4)

    plt.xlim((0, max))
    # plt.xlim((0, 10000))

    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xticks([0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000],
               [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
    # plt.xticks([0, 10000, 20000, 30000, 40000, 50000],
    #            [0, 20, 40, 60, 80, 100])
    # plt.xticks([0, 2000, 4000, 6000, 8000, 10000],
    #            [0, 4, 8, 12, 16, 20])

    # plt.xlabel("times(ms)", fontsize=13)
    plt.xlabel("time", fontsize=14)
    plt.ylabel("concentration", fontsize=14)
    date_str = time.strftime('%Y-%m-%d-%H-%M-%S', time.localtime())
    # plt.savefig("../figure/line/" + date_str + "_" + type + parameter + ".png")
    plt.savefig("a.png")
    plt.show()


if __name__ == '__main__':
    # plot(label_list, "fn", filename, title, paths)
    plot(label_list, "gn", filename, title, paths)
