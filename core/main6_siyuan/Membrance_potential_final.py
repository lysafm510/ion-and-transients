# 蟹不肉，这里是大蛋儿的代码空间
# 钙火花的信息将由大蛋儿为你展现
# skr，skr~~~

import numpy as np
import math
import os

global Ca_cyto, Ca_fsr, K_ryr, K_fsr, BCSQ, KCSQ, Dca, Df, Kf1, Kf2, Ft  # 对应文档中的 2.2新方程中的6 中常数，其中Kf1是k+，Kf2是k-
global DT  # DT为△t
global Fn, Gn, NEWFn, NEWGn
global Determinant, Area, control_area, relevant_triangle
global Cm, V_in, V_out, T, R, F, Vj


def not_empty(s):
    return s and s.strip()


def Grid():  # 网格的信息
    global NP, NE, PMAX, EMAX, NVEX
    global NPOCH, NOD, NOE, PX, PY

    NP = 0  # 点的个数
    NE = 0  # 三角形的个数
    PMAX = 50000  # 点数
    EMAX = 100000  # 三角形单元数
    NVEX = 3  # 三角形顶点数

    NPOCH = np.empty(PMAX, int)  # 每个节点的属性
    NOD = np.empty([NVEX, EMAX], int)  # 这些节点组成的三角形单元的信息
    NOE = np.empty([NVEX, EMAX], int)  # 这个是每一个三角形单元的三个邻居单元的编号
    PX = np.empty(PMAX, float)  # 点的横坐标
    PY = np.empty(PMAX, float)  # 点的纵坐标

    gridt_file = open('../main5_beibei/GRID/gridt.dat', 'r')  # gridt.dat，这个是空间离散化的每一个离散的节点的坐标，第一行是节点的数量
    gridt_count = 0
    for line in gridt_file.readlines():
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        PX[gridt_count] = current_line[0]
        PY[gridt_count] = current_line[1]
        gridt_count = gridt_count + 1
    NP = gridt_count  # NP = 6413
    gridt_file.close()

    nod_file = open('../main5_beibei/GRID/nod.dat',
                    'r')  # nod.dat，这个是gridt.dat中这些节点组成的三角形单元的信息，第一行是三角形单元的数量，后面是每一个三角元对应的三个节点的编号；
    nod_count = 0
    for line in nod_file.readlines():
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        NOD[0, nod_count] = int(current_line[0]) - 1
        NOD[1, nod_count] = int(current_line[1]) - 1
        NOD[2, nod_count] = int(current_line[2]) - 1
        nod_count = nod_count + 1
    NE = nod_count  # NE = 12444
    nod_file.close()

    noe_file = open('../main5_beibei/GRID/noe.dat',
                    'r')  # noe.dat，这个是每一个三角形单元的三个邻居单元的编号（所谓邻居指的是和这个单元有一条边重合的三角形单元），编号0表示这条边是区域边界，所以没有邻居
    noe_count = 0
    for line in noe_file.readlines():
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        NOE[0, noe_count] = current_line[0]
        NOE[1, noe_count] = current_line[1]
        NOE[2, noe_count] = current_line[2]
        noe_count = noe_count + 1
    noe_file.close()

    npoch_file = open('../main5_beibei/GRID/npoch.dat',
                      'r')  # npoch.dat为每个节点的属性：0代表内部节点，2代表入流边界或远场的节点，1代表壁面的节点，4代表出流边界节点
    npoch_count = 0
    for line in npoch_file.readlines():
        current_line = list(filter(not_empty, line.strip("\n").split(" ")))
        NPOCH[npoch_count] = current_line[0]
        npoch_count = npoch_count + 1
    npoch_file.close()
    print("执行网格代码完成！")


def Initialization():  # 参数初始化
    global Ca_cyto, Ca_fsr, K_ryr, K_fsr, BCSQ, KCSQ, Dca, Df, Kf1, Kf2, Ft
    global Fn, Gn, NEWFn, NEWGn, midfn, midgn
    global NP, NE
    global Determinant, Area, control_area, relevant_triangle, Total_area
    global Allnandl, Allabc, nmax
    global Inner, Wall, Inflow, Outflow
    global DT, RELEASE_TIMES
    global Cm, V_in, V_out, T, R, F, Vj
    Inner = 0  # 0代表内部节点，2代表入流边界或远场的节点，1代表壁面的节点，4代表出流边界节点
    Wall = 1
    Inflow = 2
    Outflow = 4
    Ca_cyto = 0.0001
    Ca_fsr = 1.0
    K_ryr = 0.25 * 6.5 * 10 ** 7
    # K_fsr = 0.785 * 10 ** 6
    K_fsr = 0.7854 * 10 ** 6
    BCSQ = 14.0
    KCSQ = 0.63
    Dca = 3.5 * 10 ** 8
    Df = 2.0 * 10 ** 7
    # Kf1 = 1000
    # Kf2 = 400
    Kf1 = 48800
    Kf2 = 19520
    Ft = 0.1

    # 关闭RyR通道用到,*****************************************************这里怎么取的
    DT = 2 * 10 ** -6  # dt
    RELEASE_TIMES = 2 * 10 ** -2  # 0.02s,这个就是ryr通道开放的时间，后80毫秒是恢复的

    Fn = np.empty(PMAX, float)
    # Gn = np.empty(PMAX, float)
    Gn = np.ones(PMAX, float) * (1 / 14)  # gn初始值 1/14  固定
    NEWFn = np.empty(PMAX, float)
    NEWGn = np.empty(PMAX, float)
    midfn = np.empty(NP, float)
    midgn = np.empty(NP, float)
    Determinant = np.zeros(NE, float)  # 计算三角形面积用到的行列式D
    Area = np.zeros(NE, float)  # 三角形的面积
    Total_area = 0.0  # 整个区域的面积，在后边求fn和gn的平均值时要用
    for i in range(0, NP):
        Fn[i] = 1.0
        # Gn[i] = 0.0909

    matrix = np.zeros((3, 3))  # 计算三角形面积用到的行列式D
    control_area = np.zeros(NP, float)  # 存储每个点的控制面积
    count = np.zeros(NP, int)
    nmax = 0
    for i in range(0, NE):  # 计算每个三角形的面积
        count[NOD[0, i]] = count[NOD[0, i]] + 1
        count[NOD[1, i]] = count[NOD[1, i]] + 1
        count[NOD[2, i]] = count[NOD[2, i]] + 1
        # matrix = np.zeros((3, 3))
        for j in range(0, 3):
            matrix[j, 0] = 1.0
            matrix[j, 1] = PX[NOD[j, i]]
            matrix[j, 2] = PY[NOD[j, i]]
        Determinant[i] = np.linalg.det(matrix)
        Area[i] = 1 / 2.0 * Determinant[i]  # np.linalg.det(A)求行列式的值,文档中面积A=1/2 * D，其中D为  Determinant
        # Area[i] = abs(1 / 2.0 * np.linalg.det(Determinant))  # 求面积需要用到绝对值吗？
        Total_area = Total_area + Area[i]
        for m in range(0, 3):
            control_area[NOD[m, i]] = control_area[NOD[m, i]] + Area[i]

    nmax = int(max(count))
    print("三角形面积计算完成！")

    Relevant_Tri()  # 初始化每个点相邻的三角形信息

    Allnandl = np.empty([NE, 3, 3], float)  # 存储每个点对应的nix, niy, Li
    Allabc = np.empty([NE, 3, 3], float)  # 存储每个点对应的a,b,c
    for i in range(0, NE):  # 求每个点为中心点时的nix, niy, Li和a,b,c
        for j in range(0, 3):
            p2 = NOD[(j + 1) % 3, i]
            p3 = NOD[(j + 2) % 3, i]
            Allabc[i, j, 0] = (PX[p2] * PY[p3] - PX[p3] * PY[p2]) / Determinant[i]  # 每个点的a
            Allabc[i, j, 1] = (PY[p2] - PY[p3]) / Determinant[i]  # 每个点的b
            Allabc[i, j, 2] = (PX[p3] - PX[p2]) / Determinant[i]  # 每个点的c

            L = math.sqrt((PX[p2] - PX[p3]) ** 2 + (PY[p2] - PY[p3]) ** 2)
            Allnandl[i, j, 0] = L  # 求每个点对应的L
            # Allnandl[i, j, 1] = (PY[p2] - PY[p3]) / L   # 求每个点对应的nx
            # Allnandl[i, j, 2] = (PX[p3] - PX[p2]) / L   # 求每个点对应的ny
            Allnandl[i, j, 1] = (PY[p3] - PY[p2]) / L  # 求每个点对应的nx
            Allnandl[i, j, 2] = (PX[p2] - PX[p3]) / L  # 求每个点对应的ny

    Cm = 2 * (Total_area * (10 ** -4))
    V_in = 0.0
    V_out = 0.0
    T = 310  # 绝对温度，单位K
    R = 8.314  # 气体常数，单位J/(mol*K)
    F = 96485.33289  # 法拉第常数，96485.33289±0.00059C/mol
    Vj = -((R * T) / (2 * F)) * math.log(Ca_fsr / Ca_cyto)  # 初始能斯特电势
    print("执行参数初始化代码完成！")


def Relevant_Tri():
    print("开始执行Relevant_Tri")
    global NP, NE, NOD, nmax
    global relevant_triangle

    if not os.path.exists("TRIANGLE_NUMBER1"):  # 将每个点相邻的三角形保存到本地本地文件中
        os.mkdir("TRIANGLE_NUMBER1")
    # relevant_triangle = np.empty([NP, NE, 1], int)  # 存储包含每个点的三角形以及是改点是三角形中的具体哪个点其中的第三维 1 保存是某个三角形中的具体的点
    relevant_triangle = np.full([NP, nmax, 2], -1)  # 存储包含每个点的三角形以及是改点是三角形中的具体哪个点其中的第三维 1 保存是某个三角形中的具体的点
    # for i in range(0, NP):
    #     for t in range(0, NE):
    #         relevant_triangle[i, t, 0] = -1
    for i in range(0, NP):
        num = 0
        if not os.path.exists("TRIANGLE_NUMBER1\\POINT" + str(i + 1) + ".dat"):
            with open("TRIANGLE_NUMBER1\\POINT" + str(i + 1) + ".dat", "w") as file_object:
                for j in range(0, NE):
                    for m in range(0, 3):
                        if (i == NOD[m, j]):
                            relevant_triangle[i, num, 0] = j
                            relevant_triangle[i, num, 1] = m
                            file_object.write(str(j) + " ")
                            file_object.write(str(m) + "\n")
                    num = num + 1
            file_object.close()
        else:
            with open("TRIANGLE_NUMBER1\\POINT" + str(i + 1) + ".dat", "r") as file_object:
                for line in file_object.readlines():
                    current_line = list(filter(not_empty, line.strip("\n").split(" ")))
                    relevant_triangle[i, num, 0] = current_line[0]
                    relevant_triangle[i, num, 1] = current_line[1]
                    num = num + 1
            file_object.close()
    print("Relevant_Tri执行完成")


def sort(N, N1, N2, N3):  # 排序函数，总的来说排序这里不是很懂！！！！！
    global n2, n3, numn2, numn3
    if (N == N1):
        n2 = N2
        n3 = N3
        numn2 = 1
        numn3 = 2
    elif (N == N2):
        n2 = N3
        n3 = N1
        numn2 = 2
        numn3 = 0
    elif (N == N3):
        n2 = N1
        n3 = N2
        numn2 = 0
        numn3 = 1
    return n2, n3, numn2, numn3


def judge_boundary(triangle):
    global Inner, Wall, Inflow, Outflow, NOD, NPOCH
    global Allnandl, Allabc
    inflownum = 0
    outflownum = 0
    InL = 0.0
    OutL = 0.0
    for i in range(0, 3):
        if (NPOCH[NOD[i, triangle]] == Inflow):
            inflownum = inflownum + 1
        if (NPOCH[NOD[i, triangle]] == Outflow):
            outflownum = outflownum + 1
    if (inflownum == 2):  # 入流点
        for i in range(0, 3):
            if (NPOCH[NOD[i, triangle]] != Inflow):
                InL = Allnandl[triangle, i, 0]
    if (outflownum == 2):
        for i in range(0, 3):
            if (NPOCH[NOD[i, triangle]] != Outflow):
                OutL = Allnandl[triangle, i, 0]
    return inflownum, outflownum, InL, OutL


def iterations_fn(iterations_num):
    global Determinant, Area, control_area, relevant_triangle
    global NEWFn, midfn
    global NPOCH, NOD, NOE, PX, PY
    global Allnandl, Allabc, nmax
    global DT
    global n2, n3, numn2, numn3
    global Cm, V_in, V_out, T, R, F, Vj
    step = iterations_num
    # for i in range(0, NP):
    #     midfn[i] = 0.0

    print("vm", V_in)
    print("nernst", Vj)
    print("DCARYR", (1 - V_in * (1 / Vj)) * K_ryr)

    while (iterations_num > 0):
        for i in range(0, NP):
            f_2 = DT * (Kf2 * Gn[i] - Kf1 * Fn[i] * (Ft - Gn[i])) * control_area[i]  # 这里的fn用n时刻的值还是每次迭代的值？？？？？？
            f_5 = 0.0
            f_1sum = 0.0
            f_4sum = 0.0
            f_7sum = 0.0
            f_8sum = 0.0
            for j in range(0, nmax):
                if relevant_triangle[i, j, 0] != -1:  # 如果该点在某个三角形中则进行计算
                    point = relevant_triangle[i, j, 1]
                    triangle = relevant_triangle[i, j, 0]
                    # for m in range(0, 3):    # 取该三角形的其他两个顶点
                    #     if(point == NOD[m,j]):
                    #         # n2 = NOD[(m + 1) % 3, j]  # 排序是否可以用这个取代？？
                    #         # n3 = NOD[(m + 2) % 3, j]
                    n2, n3, numn2, numn3 = sort(i, NOD[0, triangle], NOD[1, triangle],
                                                NOD[2, triangle])  # 对这个三角形的三个点进行排序
                    if step == iterations_num:
                        # fni = (Fn[NOD[0, triangle]] + Fn[NOD[1, triangle]] + Fn[NOD[2, triangle]]) / 3     这里有错误，在求三个点的平均值时，三个点的取值也应该符合，第一次取n时刻的值，以后每次迭代取n时刻的值和上次迭代的平均值
                        # Ca_n2 = Fn[NOD[index2, j]]   # f_2i和f_3i
                        # Ca_n3 = Fn[NOD[index3, j]]
                        Ca_n1 = Fn[i]
                        Ca_n2 = Fn[n2]  # f_2i和f_3i
                        Ca_n3 = Fn[n3]

                        Ca_n2_jk = Fn[n2]  # f_2j、f_3j和f_2k、f_3k
                        Ca_n3_jk = Fn[n3]
                    else:
                        # fni = (NEWFn[NOD[0, triangle]] + NEWFn[NOD[1, triangle]] + NEWFn[NOD[2, triangle]]) / 3
                        Ca_n1 = (NEWFn[i] + Fn[i]) / 2.0
                        Ca_n2 = (NEWFn[n2] + Fn[n2]) / 2.0
                        Ca_n3 = (NEWFn[n3] + Fn[n3]) / 2.0

                        Ca_n2_jk = NEWFn[n2]  # f_2j、f_3j和f_2k、f_3k
                        Ca_n3_jk = NEWFn[n3]
                    fni = (Ca_n1 + Ca_n2 + Ca_n3) / 3.0
                    f_5 = f_5 + (1 + (BCSQ * KCSQ) / ((KCSQ + fni) ** 2)) * Area[triangle]

                    if (NPOCH[n2] == 0 or NPOCH[
                        n3] == 0):  # 如果这个三角形中其余的两个顶点不是入流、出流、壁面边界则进行1,4,6式的计算，对应文档中优化5跳过符合上述条件的三角形
                        f_1sum = f_1sum + (Allnandl[triangle, point, 1] * (
                                    Ca_n2 * Allabc[triangle, numn2, 1] + Ca_n3 * Allabc[triangle, numn3, 1]) +
                                           Allnandl[triangle, point, 2] * (
                                                       Ca_n2 * Allabc[triangle, numn2, 2] + Ca_n3 * Allabc[
                                                   triangle, numn3, 2])) * Allnandl[triangle, point, 0]
                        f_4sum = f_4sum + (Allnandl[triangle, point, 1] * Allabc[triangle, point, 1] + Allnandl[
                            triangle, point, 2] * Allabc[triangle, point, 2]) * Allnandl[triangle, point, 0]

                    innum, outnum, inl, outl = judge_boundary(triangle)  # 判断三角形中是否含有入流和出流边界
                    if (innum == 2):  # 包含入流边界
                        f_7sum = f_7sum + K_fsr * (Ca_fsr - (Ca_n2_jk + Ca_n3_jk) / 3.0) * inl
                        f_8sum = f_8sum + K_fsr * inl
                    if (outnum == 2):  # 包含出流边界
                        f_7sum = f_7sum + (1 - V_in * (1 / Vj)) * K_ryr * (Ca_cyto - (Ca_n2_jk + Ca_n3_jk) / 3.0) * outl
                        f_8sum = f_8sum + (1 - V_in * (1 / Vj)) * K_ryr * outl
            f_3 = Fn[i] * f_5
            f_1 = DT * Dca * f_1sum
            f_4 = Fn[i] * DT * Dca * f_4sum / 2.0
            f_6 = DT * Dca * f_4sum / 2.0
            f_7 = DT * f_7sum
            f_8 = DT * f_8sum / 3.0
            midfn[i] = (f_1 + f_2 + f_3 + f_4 + f_7) / (f_5 - f_6 + f_8)
        for j in range(0, NP):
            NEWFn[j] = midfn[j]
        iterations_num = iterations_num - 1
    # for I in range(0, NP):
    #     Fn[I] = NEWFn[I]


def iterations_gn(iterations_num):
    global Determinant, Area, control_area, relevant_triangle
    global NEWGn, midgn
    global NPOCH, NOD, NOE, PX, PY
    global Allnandl, Allabc
    global DT
    global n2, n3, numn2, numn3
    step = iterations_num
    # for i in range(0, NP):
    #     NEWGn[i] = 0.0
    while (iterations_num > 0):
        for i in range(0, NP):
            g_1sum = 0.0
            g_2 = DT * (Kf1 * Fn[i] * (Ft - Gn[i]) - Kf2 * Gn[i]) * control_area[i]
            g_3 = Gn[i] * control_area[i]
            g_4sum = 0.0
            g_5 = control_area[i]

            for j in range(0, nmax):
                if relevant_triangle[i, j, 0] != -1:  # 如果该点在某个三角形中则进行计算
                    point = relevant_triangle[i, j, 1]
                    triangle = relevant_triangle[i, j, 0]
                    # for m in range(0, 3):    # 取该三角形的其他两个顶点
                    #     if(point == NOD[m,j]):
                    n2, n3, numn2, numn3 = sort(i, NOD[0, triangle], NOD[1, triangle],
                                                NOD[2, triangle])  # 对这个三角形的三个点进行排序
                    if step == iterations_num:
                        Gn_n2 = Gn[n2]  # g_2i和g_3i
                        Gn_n3 = Gn[n3]
                    else:
                        Gn_n2 = (NEWGn[n2] + Gn[n2]) / 2.0
                        Gn_n3 = (NEWGn[n3] + Gn[n3]) / 2.0
                    if (NPOCH[n2] == 0 or NPOCH[
                        n3] == 0):  # 如果这个三角形中其余的两个顶点不是入流、出流、壁面边界则进行1,4,6式的计算，对应文档中优化5跳过符合上述条件的三角形
                        g_1sum = g_1sum + (Allnandl[triangle, point, 1] * (
                                    Gn_n2 * Allabc[triangle, numn2, 1] + Gn_n3 * Allabc[triangle, numn3, 1]) +
                                           Allnandl[triangle, point, 2] * (
                                                       Gn_n2 * Allabc[triangle, numn2, 2] + Gn_n3 * Allabc[
                                                   triangle, numn3, 2])) * Allnandl[triangle, point, 0]
                        g_4sum = g_4sum + (Allnandl[triangle, point, 1] * Allabc[triangle, point, 1] + Allnandl[
                            triangle, point, 2] * Allabc[triangle, point, 2]) * Allnandl[triangle, point, 0]
            g_1 = DT * Df * g_1sum
            g_4 = DT * Gn[i] * Df * g_4sum / 2.0
            g_6 = DT * Df * g_4sum / 2.0
            midgn[i] = (g_1 + g_2 + g_3 + g_4) / (g_5 - g_6)
        for j in range(0, NP):
            NEWGn[j] = midgn[j]
        iterations_num = iterations_num - 1
    # for I in range(0, NP):
    #     Gn[I] = NEWGn[I]


def BeginLoop():
    global NP, Fn, Gn, NEWFn, midfn, midgn
    global average_fn, average_gn
    global ICARYR, ICAFSR, K_ryr
    global DT, RELEASE_TIMES
    global LENGTH_out_sum, fn_out_sum
    global Cm, V_in, V_out, T, R, F, Vj
    Currentstep = 1  # 当前步
    Loopnum = 50000  # 循环执行的次数

    iterations_num = 10  # 迭代次数
    savestep = 1  # 每几步保存一次4··

    path_savefn = "DATA_membrance_potential\\final1" + "\\Fn"  # fn保存路径
    path_savegn = "DATA_membrance_potential\\final1" + "\\Gn"  # gn保存路径
    if not os.path.exists(path_savefn):
        os.makedirs(path_savefn)
    if not os.path.exists(path_savegn):
        os.makedirs(path_savegn)

    if Currentstep == 1:  # 如果是从第一步开始
        Fn_file = path_savefn + "\\Fn00000000.dat"
        Fn_save = open(Fn_file, "w")
        Gn_file = path_savegn + "\\Gn00000000.dat"
        Gn_save = open(Gn_file, "w")
        for i in range(0, NP):
            Fn_save.write(str(Fn[i]) + "\n")
            Gn_save.write(str(Gn[i]) + "\n")
        average_Fn()  # 求平均钙离子浓度
        CA_CURRENT()  # 求Ca电流
        Vj = -((R * T) / (2 * F)) * math.log((CA_OUT / LENGTH_out_sum) / Ca_cyto)
        V_in = (- (1 / Cm) * ICARYR * DT) * (10.0 ** -3) + V_in
        Fn_save.write(str(average_fn) + "\n")  # NP + 1， 纪录当前的平均钙离子浓度
        # Fn_save.write(str() + "\n")  # NP + 2, 纪录当前步数
        Fn_save.write(str(ICARYR) + " " + str(ICAFSR) + "\n")  # NP + 2, 纪录钙离子电流
        Fn_save.write(str(Currentstep) + "\n")  # NP + 3, 纪录当前步数
        average_Gn()  # 求平均荧光钙浓度
        Gn_save.write(str(average_gn) + "\n")  # NP + 1， 纪录当前的平均荧光钙浓度
        Gn_save.write(str(Currentstep) + "\n")  # NP + 2, 纪录当前步数
        print("数据写入  " + Fn_file)
        print("数据写入  " + Gn_file)
        Fn_save.close()
        Gn_save.close()
    else:  # 如果不是从第一步开始
        step = str(Currentstep - 1).zfill(8)
        read_fn = path_savefn + "\\Fn" + step + ".dat"
        read_gn = path_savegn + "\\Gn" + step + ".dat"
        read_fn_file = open(read_fn, "r")
        read_gn_file = open(read_gn, "r")
        fnnum = 0
        gnnum = 0
        for line in read_fn_file.readlines():
            if fnnum < NP:
                Fn[fnnum] = float(line)  # 存在Fn中
            else:
                if fnnum == NP:
                    average_fn = float(line)
                elif fnnum == NP + 1:
                    tmp = line.split()
                    ICARYR = float(tmp[0])
                    ICAFSR = float(tmp[1])
                else:
                    STEP1 = int(line)
                    print("读取原有的Save", STEP1)
            fnnum = fnnum + 1
        read_fn_file.close()

        for line in read_gn_file.readlines():
            if gnnum < NP:
                Gn[gnnum] = float(line)
            else:
                if gnnum == NP:
                    average_gn = float(line)
                else:
                    STEP2 = int(line)
                    print("读取原有的Gn", STEP2)
            gnnum = gnnum + 1
        read_gn_file.close()

    RELEASE_STEP = int((RELEASE_TIMES + (10 ** -6)) / DT)  # 释放时间的迭代次数
    if Currentstep >= RELEASE_STEP:
        K_ryr = 0
    # 释放的点就停止释放

    for i in range(Currentstep, Loopnum + 1):
        for j in range(0, NP):
            NEWFn[j] = 0.0
            NEWGn[j] = 0.0
            midfn[j] = 0.0
            midgn[j] = 0.0
        iterations_fn(iterations_num)
        iterations_gn(iterations_num)
        for I in range(0, NP):
            Fn[I] = NEWFn[I]
            Gn[I] = NEWGn[I]
        if (i % savestep == 0) or (i == Loopnum):
            average_Fn()  # 求平均钙离子浓度
            average_Gn()  # 求平均荧光钙浓度
            CA_CURRENT()  # 还要求那个什么Ca电流
            Vj = -((R * T) / (2 * F)) * math.log((CA_OUT / LENGTH_out_sum) / Ca_cyto)
            V_in = (- (1 / Cm) * ICARYR * DT) * (10.0 ** -3) + V_in
            stepname = str(i).zfill(8)
            fn_file = path_savefn + "\\Fn" + stepname + ".dat"
            gn_file = path_savegn + "\\Gn" + stepname + ".dat"
            savefn = open(fn_file, "w")
            savegn = open(gn_file, "w")
            for j in range(0, NP):
                savefn.write(str(Fn[j]) + "\n")
                savegn.write(str(Gn[j]) + "\n")
            savefn.write(str(average_fn) + "\n")  # NP + 1， 纪录当前的平均钙离子浓度
            # Fn_save.write(str(ICARYR) + " " + str(ICAFSR) + "\n")  # NP + 2, 纪录钙离子电流
            savefn.write(str(ICARYR) + " " + str(ICAFSR) + "\n")  # NP + 2, 纪录钙离子电流
            savefn.write(str(i) + "\n")  # NP + 3, 纪录当前步数
            savegn.write(str(average_gn) + "\n")  # NP + 1， 纪录当前的平均荧光钙浓度
            savegn.write(str(i) + "\n")  # NP + 2, 纪录当前步数
            print("数据写入 ", fn_file)
            print("数据写入 ", gn_file)
            savefn.close()
            savegn.close()

        if i == RELEASE_STEP:
            K_ryr = 0
            #
            # if J >= 50:
            #     DCARYR = 0

            # print("DCARYR ", DCARYR)
            # for i in range(0, NP):
            #     NEWCCA[i] = 0.
            #     NEWF[i] = 0.

        print("执行成功")


def average_Fn():
    """
    求平均钙离子浓度
    """
    global total_fn, average_fn
    global Fn, control_area, Total_area
    total_fn = 0.0
    for i in range(0, NP):
        total_fn = total_fn + Fn[i] * control_area[i]  # 每一点Ca的浓度X每一个点对应的三角形单元面积=每一点钙的总面积
    total_fn = total_fn / 3
    if Total_area > 0.00000000000001:
        average_fn = total_fn / Total_area  # JSR平均浓度


def average_Gn():
    """
    求平均荧光钙浓度
    """
    global total_gn, average_gn
    global Gn, control_area, Total_area
    total_gn = 0.0
    for i in range(0, NP):
        total_gn = total_gn + Gn[i] * control_area[i]  # 每一点Ca的浓度X每一个点对应的三角形单元面积=每一点钙的总面积
    total_gn = total_gn / 3
    if Total_area > 0.00000000000001:
        average_gn = total_gn / Total_area  # JSR平均浓度


def CA_CURRENT():
    """
    求钙离子电流
    """
    global ICARYR, ICAFSR, average_fn, average_gn, LENGTH
    global Inflow, Outflow
    global K_ryr, K_fsr
    global LENGTH_out_sum, fn_out_sum, CA_OUT
    H_JSR = 30.0  # JSR的高度
    UNITEC = 1.610217733  # 一个电荷的常数
    MOLNUM = 6.0221367  # 摩尔常数，单位
    ICARYR = 0.0
    ICAFSR = 0.0
    CA_IN = 0.0
    CA_OUT = 0.0
    LENGTH = 0.0
    LENGTH_out_sum = 0.0
    # fn_out_sum = 0.0
    for I in range(0, NE):
        N1 = NOD[0, I]
        N2 = NOD[1, I]
        N3 = NOD[2, I]
        #  OUT CURRENT
        if (NPOCH[N1] == Outflow) and (NPOCH[N2] == Outflow):
            LENGTH = math.sqrt((PX[N1] - PX[N2]) ** 2 + (PY[N1] - PY[N2]) ** 2)
            CA_OUT = CA_OUT + LENGTH * (Fn[N1] + Fn[N2]) / 2
            LENGTH_out_sum = LENGTH_out_sum + LENGTH
            # fn_out_sum = fn_out_sum + Fn[N1] + Fn[N2]
        if (NPOCH[N3] == Outflow) and (NPOCH[N2] == Outflow):
            LENGTH = math.sqrt((PX[N3] - PX[N2]) ** 2 + (PY[N3] - PY[N2]) ** 2)
            CA_OUT = CA_OUT + LENGTH * (Fn[N3] + Fn[N2]) / 2
            LENGTH_out_sum = LENGTH_out_sum + LENGTH
            # fn_out_sum = fn_out_sum + Fn[N3] + Fn[N2]
        if (NPOCH[N1] == Outflow) and (NPOCH[N3] == Outflow):
            LENGTH = math.sqrt((PX[N1] - PX[N3]) ** 2 + (PY[N1] - PY[N3]) ** 2)
            CA_OUT = CA_OUT + LENGTH * (Fn[N1] + Fn[N3]) / 2
            LENGTH_out_sum = LENGTH_out_sum + LENGTH
            # fn_out_sum = fn_out_sum + Fn[N1] + Fn[N3]
        #  IN CURRENT
        if (NPOCH[N1] == Inflow) and (NPOCH[N2] == Inflow):
            LENGTH = math.sqrt((PX[N1] - PX[N2]) ** 2 + (PY[N1] - PY[N2]) ** 2)
            CA_IN = CA_IN + LENGTH * (Ca_fsr - (Fn[N1] + Fn[N2]) / 2)
        if (NPOCH[N3] == Inflow) and (NPOCH[N2] == Inflow):
            LENGTH = math.sqrt((PX[N3] - PX[N2]) ** 2 + (PY[N3] - PY[N2]) ** 2)
            CA_IN = CA_IN + LENGTH * (Ca_fsr - (Fn[N3] + Fn[N2]) / 2)
        if (NPOCH[N1] == Inflow) and (NPOCH[N3] == Inflow):
            LENGTH = math.sqrt((PX[N1] - PX[N3]) ** 2 + (PY[N1] - PY[N3]) ** 2)
            CA_IN = CA_IN + LENGTH * (Ca_fsr - (Fn[N1] + Fn[N3]) / 2)

    CA_OUT = CA_OUT * (1 - V_in * (1 / Vj)) * K_ryr * H_JSR
    CA_IN = CA_IN * K_fsr * H_JSR
    ICARYR = CA_OUT * UNITEC * 2 * MOLNUM * (10.0 ** -11)
    ICAFSR = CA_IN * UNITEC * 2 * MOLNUM * (10.0 ** -11)


def main():
    print("开始运行")
    print("运行网格函数")
    Grid()  # 网格的信息
    print("运行初始化函数")
    Initialization()  # 参数初始化
    # return
    print("循环函数")
    BeginLoop()  # 开始循环


"""
保佑代码没bug
"""
main()  # 保佑代码没bug

"""
保佑代码没bug
"""

# 第一次运行到612
# 第二次运行到1926
# 第三次运行到2181
# 第四次运行到4338,   这一次为了减少时间复杂度，合并了几个for循环
# 第五次运行到4490
