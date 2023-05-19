import os

import pandas as pd

DATA_PATH = "test15_20220707_center4_new_area_no"
types_fn = ['avg_ca_jsr', 'c_ca_store', 'i_ca_ryr', 'i_ca_fsr', 'vm', 'nernst', 'd_ca_ryr_jsr']
types_gn = ['avg_gn_jsr']


def extract(types, pattern):
    for filename in os.listdir("../data/" + DATA_PATH + "/blink/" + pattern):
        df = pd.read_csv("../data/" + DATA_PATH + "/blink/" + pattern + "/" + filename, dtype=str, header=None,
                         index_col=0)
        for type in types:
            with open("../data/" + DATA_PATH + "/" + type + ".csv", "a") as f:
                f.writelines(str(float(df.loc[type])) + '\n')


if __name__ == '__main__':
    blink_fn_path = "../../data/" + DATA_PATH + "/blink/fn"
    blink_gn_path = "../../data/" + DATA_PATH + "/blink/gn"
    # extract(types_fn, "fn")
    extract(types_gn, "gn")
