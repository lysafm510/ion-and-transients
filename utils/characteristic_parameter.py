import pandas as pd

path = "test15_20220707_center4_new_area_store_5"

data_frame = pd.read_csv("../data/" + path + "/avg_gn_jsr.csv", header=None)
data = data_frame[0].tolist()


def cal_amplitude():
    # ΔF/F0
    print("最大值：%s" % max(data))
    print("最小值：%s" % min(data))
    print("最小值步数：%s" % (data.index(min(data)) * 100))

    amplitude = (max(data) - min(data)) / max(data)
    print("振幅ΔF/F0：%s" % amplitude)
    print("******************************************* ")
    print()


def cal_t50():
    dt = 2 * 10 ** -6
    peak_time = (data.index(min(data))) * 100 * dt
    half_concentration = (max(data) - min(data)) / 2 + min(data)
    print("tRise：%s" % peak_time)
    print("峰值一半时浓度：%s" % half_concentration)

    index = 100
    while index < len(data) - 1 and data[index] < half_concentration:
        index = index + 1
    if index >= len(data) - 1:
        print("最后一步浓度：%s" % data[index])
        print("暂未回升到一半")
    else:
        left_concentration = data[index - 1]
        right_concentration = data[index]
        half_time = index * 100 * dt - (right_concentration - half_concentration) * dt / (
                right_concentration - left_concentration)
        print("t50：%s" % ((half_time - peak_time) * 1000))
    print("******************************************* ")
    print()


def cal_FDHM():
    dt = 2 * 10 ** -6
    half_concentration = (max(data) - min(data)) / 2 + min(data)

    index = 0
    while index < len(data) - 1 and data[index] > half_concentration:
        index = index + 1
    left_concentration = data[index - 1]
    right_concentration = data[index]
    half_time1 = index * 100 * dt - (right_concentration - half_concentration) * dt / (
            right_concentration - left_concentration)
    print("下降到一半时间：%s" % half_time1)

    index = 100
    while index < len(data) - 1 and data[index] < half_concentration:
        index = index + 1
    if index >= len(data) - 1:
        print("暂未回升到一半")
    else:
        left_concentration = data[index - 1]
        right_concentration = data[index]
        half_time2 = index * 100 * dt - (right_concentration - half_concentration) * dt / (
                right_concentration - left_concentration)
        print("回升到一半时间：%s" % half_time2)
        print("FDHM：%s" % ((half_time2 - half_time1) * 1000))


if __name__ == '__main__':
    if len(data) < 100:
        print("暂未下降到最低点")
    else:
        cal_amplitude()
        cal_t50()
        cal_FDHM()
