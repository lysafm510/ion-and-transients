import random
from math import sqrt, pi, exp, log

import numpy as np
import pandas as pd

custom_arr = [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62,
              64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]
custom_arr = np.array(custom_arr)
a = 200 - custom_arr
b = np.append(custom_arr,a)
b = np.sort(b)
# print(a)
print(b.tolist())
