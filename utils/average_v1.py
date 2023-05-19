import numpy as np
import linecache
import os

path = "test15_20220707_center4_new_area_no"


def get_avg_line():
    data = np.loadtxt("../data/" + path + "/blink/fn/Fn_00000000.csv")
    # return len(data) - 4
    return len(data)


def average(line, file_name, pattern):
    with open("../data/" + path + "/" + file_name, "w") as file_object:
        for filename in os.listdir("../data/" + path + "/blink/" + pattern):
            data = linecache.getline("../data/" + path + "/blink/" + pattern + "/" + filename, line)
            file_object.write(data)


if __name__ == '__main__':
    avg_line = get_avg_line()

    # average(avg_line - 2, "d_ca_ryr.csv", "fn")
    # average(avg_line - 3, "nernst.csv", "fn")
    # average(avg_line - 4, "vm.csv", "fn")
    # average(avg_line - 5, "i_ca_fsr.csv", "fn")
    # average(avg_line - 6, "i_ca_ryr.csv", "fn")
    # average(avg_line - 7, "avg_fn_jsr.csv", "fn")
    #
    #
    # average(avg_line - 7, "avg_gn_jsr.csv", "gn")



    average(avg_line - 2, "d_ca_ryr.csv", "fn")
    average(avg_line - 3, "nernst.csv", "fn")
    average(avg_line - 4, "vm.csv", "fn")
    average(avg_line - 5, "i_ca_fsr.csv", "fn")
    average(avg_line - 6, "i_ca_ryr.csv", "fn")
    average(avg_line - 7, "avg_fn_jsr.csv", "fn")
    average(avg_line - 8, "c_ca_store.csv", "fn")

    average(avg_line - 8, "avg_gn_jsr.csv", "gn")
