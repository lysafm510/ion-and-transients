import csv
import os
from math import fmod, log

import numpy as np
import pandas as pd

from blink import average_concentration, ca_current, blink_fn_equation, blink_gn_equation, cca_ryr_inside_store
from constant import SAVE_PATH, CCAJSR, CCAF, START_STEP, RELEASE_TIME, DT, END_TIME, SAVE_INTERVAL, DCARYR, CCAMYO, F, \
    Z, T, R, C
from grid_info import NP
from blink_coeff_matrix import TOTAL_AREA

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
    vm = 0

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
            writer.writerow([vm])
            writer.writerow([i_ca_ryr])
            writer.writerow([i_ca_fsr])
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
        fn = np.array(f_data.iloc[:len(f_data) - 4, 0]).astype('float64')
        print("数据读取Fn_" + str(START_STEP - 1).zfill(8) + ".csv...")

        g_data = pd.read_csv(blink_gn_path + "/Gn_" + str(START_STEP - 1).zfill(8) + ".csv", dtype=str, header=None)
        gn = np.array(g_data.iloc[:len(g_data) - 2, 0]).astype('float64')
        print("数据读取Gn_" + str(START_STEP - 1).zfill(8) + ".csv...")

        current_step = START_STEP - 1
        times = float(f_data.iloc[(len(f_data) - 1):, :].values[0])
        print("数据读取Spark_" + str(START_STEP - 1).zfill(8) + ".csv...")

        vm = float(f_data.iloc[(len(f_data) - 4):, :].values[0])
        i_ca_ryr = float(f_data.iloc[(len(f_data) - 3):, :].values[0])

    if times >= (RELEASE_TIME - DT / 2):
        d_ca_ryr_jsr = 0
        print("关闭RyR通道")

    while times <= END_TIME:
        times = times + DT
        current_step = current_step + 1

        coeff = - R * T / (Z * F)
        if not d_ca_ryr_jsr == 0:
            # 计算当前膜电位
            vm = (vm - i_ca_ryr * DT / (C * TOTAL_AREA)) * (10.0 ** 11)
            # 计算平衡电位
            print("vm", vm)
            c_ca_store = cca_ryr_inside_store(fn)
            print(c_ca_store)
            nernst = coeff * log(c_ca_store / CCAMYO) * (10 ** 3)
            print("nernst", nernst)
            # 计算DCARYR系数
            f = 1 - vm / nernst
            print("f", vm / nernst)
            print("1-f", f)
        else:
            f = 0

        # 计算blink
        new_fn = blink_fn_equation(fn, gn, CCAMYO, f * d_ca_ryr_jsr, 10)
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
                writer.writerow([vm])
                writer.writerow([i_ca_ryr])
                writer.writerow([i_ca_fsr])
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