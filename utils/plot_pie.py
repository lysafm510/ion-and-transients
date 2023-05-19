import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# 坐标位置
grid_path = "beibei"
# 浓度值位置
value_type = "gn"
value_path = "data/test10_测试gn上升_带电荷力/blink/gn/Gn_00000100.csv"
# 保存位置
save_path = "test"

# *****************************************************************************
point = np.loadtxt("../parameters/grid/blink/" + grid_path + "/gridt.dat", dtype=np.float64)
value_df = pd.read_csv("../" + value_path, header=None)
value = 0
if value_type == "fn":
    value = np.array(value_df.iloc[:len(value_df) - 4, 0]).astype('float64')
elif value_type == "gn":
    value = np.array(value_df.iloc[:len(value_df) - 3, 0]).astype('float64')
x, y = np.mgrid[-300:300:1000j, -300:300:1000j]  # 生成网格

grid = griddata(point, value, (x, y), method='linear')  # 对z插值

plt.contourf(x, y, grid, cmap=plt.cm.summer)
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("../figure/pie/" + save_path + ".png")
plt.show()
