import csv
import os
from math import fmod

import numpy as np
import pandas as pd

from blink import average_concentration, ca_current, blink_fn_equation, blink_gn_equation
from constant import SAVE_PATH, CCAJSR, CCAF, START_STEP, RELEASE_TIME, DT, END_TIME, SAVE_INTERVAL, DCARYR, CCAMYO
from grid_info import NP

'''
耦联子
'''


def do_loop_diffusion():
    # 保存路径
    blink_fn_path = "../../data/" + SAVE_PATH + "/blink/fn"
    blink_gn_path = "../../data/" + SAVE_PATH + "/blink/gn"

    # 临时变量，0.02s后置0
    d_ca_ryr_jsr = DCARYR

    # blink赋初值
    fn = np.full(NP, CCAJSR)  # fn,肌质网每一点钙的浓度,为肌质网每一点的初始钙浓度赋值为1
    gn = np.full(NP, CCAF)  # gn,初值1/14

    if START_STEP == 0:  # 开始前，先存储第0步初值
        times = 0.0  # 当前时间
        current_step = 0  # 当前步数

        # 写初值
        blink_fn_folder = os.path.exists(blink_fn_path)
        if not blink_fn_folder:
            os.makedirs(blink_fn_path)
        with open(blink_fn_path + "/Fn_00000000.csv", "w", newline='') as file:
            writer = csv.writer(file)
            for i in range(0, NP):
                writer.writerow([fn[i]])
            avg_ca_jsr = average_concentration(fn)
            i_ca_ryr, i_ca_fsr = ca_current(fn, d_ca_ryr_jsr)
            writer.writerow([avg_ca_jsr])
            writer.writerow([i_ca_ryr])
            writer.writerow([i_ca_fsr])
            writer.writerow([current_step])
            writer.writerow([times])
        file.close()
        print("数据写入 Fn_00000000.csv...")

        blink_gn_folder = os.path.exists(blink_gn_path)
        if not blink_gn_folder:
            os.makedirs(blink_gn_path)
        with open(blink_gn_path + "/Gn_00000000.csv", "w", newline='') as file:
            writer = csv.writer(file)
            for i in range(0, NP):
                writer.writerow([gn[i]])
            avg_gn_jsr = average_concentration(gn)
            writer.writerow([avg_gn_jsr])
            writer.writerow([times])
        file.close()
        print("数据写入Gn_00000000.csv...")

    else:  # 不是第一步开始，先读取先前一步保存的结果
        f_data = pd.read_csv(blink_fn_path + "/Fn_" + str(START_STEP - 1).zfill(8) + ".csv", dtype=str, header=None)
        fn = np.array(f_data.iloc[:len(f_data) - 5, 0]).astype('float64')
        current_step = int(f_data.iloc[(len(f_data) - 2), :].values[0])
        times = float(f_data.iloc[(len(f_data) - 1), :].values[0])
        print("数据读取Fn_" + str(START_STEP - 1).zfill(8) + ".csv...")

        g_data = pd.read_csv(blink_gn_path + "/Gn_" + str(START_STEP - 1).zfill(8) + ".csv", dtype=str, header=None)
        gn = np.array(g_data.iloc[:len(g_data) - 2, 0]).astype('float64')
        print("数据读取Gn_" + str(START_STEP - 1).zfill(8) + ".csv...")

    if times >= (RELEASE_TIME - DT / 2):
        d_ca_ryr_jsr = 0
        print("关闭RyR通道")

    while times <= END_TIME:
        times = times + DT
        current_step = current_step + 1

        # 计算blink
        new_fn = blink_fn_equation(fn, gn, CCAMYO, d_ca_ryr_jsr, 10)
        gn = blink_gn_equation(fn, gn, 10)
        for i in range(0, NP):
            fn[i] = new_fn[i]

        if (fmod(current_step, SAVE_INTERVAL) == 0) or (times == END_TIME):
            with open(blink_fn_path + "/Fn_" + str(current_step).zfill(8) + ".csv", "w", newline='') as file:
                writer = csv.writer(file)
                for i in range(0, NP):
                    writer.writerow([fn[i]])
                avg_ca_jsr = average_concentration(fn)
                i_ca_ryr, i_ca_fsr = ca_current(fn, d_ca_ryr_jsr)
                writer.writerow([avg_ca_jsr])
                writer.writerow([i_ca_ryr])
                writer.writerow([i_ca_fsr])
                writer.writerow([current_step])
                writer.writerow([times])
            file.close()
            print("数据写入Fn_" + str(current_step).zfill(8) + ".csv...")

            with open(blink_gn_path + "/Gn_" + str(current_step).zfill(8) + ".csv", "w", newline='') as file:
                writer = csv.writer(file)
                for i in range(0, NP):
                    writer.writerow([gn[i]])
                avg_gn_jsr = average_concentration(gn)
                writer.writerow([avg_gn_jsr])
                writer.writerow([times])
            file.close()
            print("数据写入Gn_" + str(current_step).zfill(8) + ".csv...")

        if abs(times - RELEASE_TIME) <= DT / 2:
            d_ca_ryr_jsr = 0
            print("关闭RyR通道")


print(SAVE_PATH)
do_loop_diffusion()
