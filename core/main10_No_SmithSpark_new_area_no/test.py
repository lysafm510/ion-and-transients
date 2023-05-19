import numpy as np
import pandas as pd

from constant import SAVE_PATH

blink_fn_path = "../../data/" + SAVE_PATH + "/blink/fn"
blink_gn_path = "../../data/" + SAVE_PATH + "/blink/gn"
START_STEP = 2
GRID_PATH = "25RyR_radius_0.8_center4_0501"
grid = np.loadtxt('../../parameters/grid/blink/' + GRID_PATH + '/gridt.dat', dtype=np.float64)
NP = len(grid)

f_data = pd.read_csv(blink_fn_path + "/Fn_" + str(START_STEP - 1).zfill(8) + ".csv", dtype=str, header=None, index_col=0)
fn = np.array(f_data.iloc[0:NP, 0]).astype('float64')
# avg_ca_jsr1 = float(f_data.iloc[(len(f_data) - 7), :].values[1])

c = float(f_data.loc['avg_ca_jsr'])
i_ca_ryr = float(f_data.loc['i_ca_ryr'])
i_ca_fsr = float(f_data.loc['i_ca_fsr'])
vm = float(f_data.loc['vm'])
nernst = float(f_data.loc['nernst'])


# print(avg_ca_jsr1)
print(f_data)

print(i_ca_ryr)
print(i_ca_fsr)
print(vm)
print(nernst)