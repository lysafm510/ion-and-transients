import numpy as np
import pandas as pd

f_data = pd.read_csv("Fn_00001000.csv", dtype=str, header=None)
fn = np.array(f_data.iloc[:len(f_data) - 4, 0]).astype('float64')
times = np.array(f_data.iloc[(len(f_data) - 4), :])[0]
print(fn)
print(times)
