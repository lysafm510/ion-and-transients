import numpy as np

from constant import GRID_PATH

"""
加载Blink网格信息
"""
# gridt.dat，这个是空间离散化的每一个离散的节点的坐标，第一行是节点的数量

grid = np.loadtxt('../../parameters/grid/blink/' + GRID_PATH + '/gridt.dat', dtype=np.float64)
# 把空间里面离散的每一点的坐标的横坐标和纵坐标分别加载到PX数组和PY数组，把结点的数量存入到NP中
PX = grid[:, 0]  # 点的横坐标
PY = grid[:, 1]  # 点的纵坐标
NP = len(grid)

# nod.dat，这个是gridt.dat中这些节点组成的三角形单元的信息，
# 第一行是三角形单元的数量，后面是每一个三角元对应的三个节点的编号
# 把每一个三角形单元所对应的三角形单元的编号分别存入到数组NOD中，即NOD是每个三角形单元所对应的节点编号的数组
nod = np.loadtxt('../../parameters/grid/blink/' + GRID_PATH + '/nod.dat', dtype=np.int64)
# NE为三角形单元的数量
NE = len(nod)
NOD = nod.T - 1

# noe.dat，这个是每一个三角形单元的三个邻居单元的编号（所谓邻居指的是和这个单元有一条边重合的三角形单元），
# 编号0表示这条边是区域边界，所以没有邻居
noe = np.loadtxt('../../parameters/grid/blink/' + GRID_PATH + '/noe.dat', dtype=np.int64)
NOE = noe.T

# npoch.dat为每个节点的属性：0代表内部节点，2代表入流边界或远场的节点，1代表壁面的节点，4代表出流边界节点
NPOCH = np.loadtxt('../../parameters/grid/blink/' + GRID_PATH + '/npoch.dat', dtype=np.int64)

# *******************************************************************************************
"""
加载SmithSpark网格信息
"""
RADIUS = np.loadtxt('../../parameters/grid/smith_spark/spark_mesh.dat', dtype=np.float64)
