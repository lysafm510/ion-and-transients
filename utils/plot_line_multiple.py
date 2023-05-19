import matplotlib.pyplot as plt
import csv
from matplotlib import font_manager

font_manager.fontManager.addfont('../fonts/time-simsun.ttf')

plt.rcParams['font.sans-serif'] = 'Times New Roman + Simsun'
plt.rcParams['axes.unicode_minus'] = False
# plt.rcParams["font.family"] = "Microsoft Yahei"
# plt.rcParams['font.sans-serif'] = ['SimHei']

path1 = "test15_20220707_center4_new_area_no"
path2 = "test15_20220707_center4_new_area_store_5"

paths = [path1, path2]

label_list = ["不考虑电荷力", "考虑电荷力"]
# 文件描述
filename = "test"
# title = "求膜电位时电流用(i_ryr-i_fsr)/2"
title = ""

save_steps_arr = [100, 100, 100]


# *****************************************************************************
def plot(label_list, type, filename, title, paths, save_steps_arr):
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
    # save_steps = 1
    for i in range(0, n):
        x = []
        y = []
        j = 0

        save_steps = save_steps_arr[i]
        with open("../data/" + paths[i] + "/avg_" + type + "_jsr.csv", 'r') as ca_csvfile1:
            plots = csv.reader(ca_csvfile1, delimiter='\t')
            for row in plots:
                x.append(j * save_steps)
                y.append(float(row[0]))
                j = j + 1
        plt.plot(x, y, linewidth=1)

        if j * save_steps > max:
            max = j * save_steps
    # plt.title(type + "   " + title, fontsize=14)
    if type == "ca":
        type1 = "[Ca]"
    else:
        type1 = "[CaF]"
    plt.title(type1, fontsize=14)
    plt.legend(loc=4, labels=label_list, fontsize=12)
    # plt.xlim((0, 50000))
    plt.xlim((0, max))
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    # plt.xticks([0, 10000, 20000, 30000, 40000, 50000], [0, 20, 40, 60, 80, 100])
    plt.xticks([0, 10000, 20000, 30000, 40000, 50000],
               [0, 20, 40, 60, 80, 100])
    # plt.xlabel("time(ms)", fontsize=13)
    plt.xlabel("time", fontsize=14)
    plt.ylabel("concentration", fontsize=14, labelpad=8.5)
    plt.grid(linestyle="--")
    # plt.savefig("../figure/test8_对冲/" + filename + "_" + type + ".png")
    plt.savefig(filename + "_" + type + ".png")
    plt.show()


if __name__ == '__main__':
    plot(label_list, "ca", filename, title, paths, save_steps_arr)
    plot(label_list, "gn", filename, title, paths, save_steps_arr)
