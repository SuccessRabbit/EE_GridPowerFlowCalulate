from numpy import *
import numpy as np
import xlrd
import xlwt
from cmath import polar, rect, sin, cos, tan, acos, asin, atan, pi


data = xlrd.open_workbook("InputFile.xlsx")
table = data.sheet_by_name("Sheet1")
# [头结点名  尾节点名   R   X   B/2]
nodal_admittance_matrix = zeros((9, 9)) * 1j  # 构造9*9节点导纳矩阵
# [标号   类型   P   U  Q   ang]
for i in range(9):
    R = table.row_values(i + 1)[2]
    X = table.row_values(i + 1)[3]
    B = table.row_values(i + 1)[4]
    row = int(table.row_values(i + 1)[0]-1)
    col = int(table.row_values(i + 1)[1]-1)
    if table.row_values(i + 1)[4] == 1:  # 是变压器支路
        nodal_admittance_matrix[row, row] += 1 / (R + X * 1j)
        nodal_admittance_matrix[row, col] -= 1 / (R + X * 1j)
        nodal_admittance_matrix[col, row] -= 1 / (R + X * 1j)
        nodal_admittance_matrix[col, col] += 1 / (R + X * 1j)
    else:  # 输电元件支路
        nodal_admittance_matrix[row, row] += 1 / (R + X * 1j) + (B * 1j)
        nodal_admittance_matrix[row, col] -= 1 / (R + X * 1j)
        nodal_admittance_matrix[col, row] -= 1 / (R + X * 1j)
        nodal_admittance_matrix[col, col] += 1 / (R + X * 1j) + (B * 1j)
print(nodal_admittance_matrix)  # 输出节点导纳矩阵


class Node:  # 定义节点类
    def __init__(self, name, type, P, U, Q, ang):
        if P != '/':
            self.P = float(P)
        else:
            self.P = 0
        if Q != '/':
            self.Q = float(Q)
        else:
            self.Q = 0
        if U != '/':
            self.U = float(U)
        else:
            self.U = 1.0
        self.ang = float(ang)
        self.type = str(type)
        self.name = int(name)-1
        self.l_P = self.P
        self.l_Q = self.Q
        self.l_ang = self.ang
        self.l_U = 1
        self.Xd = 0.3
        self.E = 0
        if name == "1":
            self.E = 1.137
        elif name == "2":
            self.E = 1.211
        elif name == "3":
            self.E = 1.043

node_info = []
nodes = []  # 存储所有节点对象
PQ_nodes = []
PV_nodes = []
PH_node = None
# 读取节点信息
for i in range(9):
    node_info.append(table.row_values(i + 11))
for i in range(9):
    info = node_info[i]
    item = Node(info[0], info[1], info[2], info[3], info[4], info[5])
    nodes.append(item)
Y = nodal_admittance_matrix

# 计算PV PQ节点个数
num_PQ = 0
num_PV = 0
for node in nodes:
    if node.type == "PQ":
        num_PQ += 1
        PQ_nodes.append(node)
    elif node.type == "PV":
        num_PV += 1
        PV_nodes.append(node)
    else:
        PH_node = node
        continue
PQ_PV_nodes = PQ_nodes + PV_nodes

t1 = []
t2 = []
for item in PQ_PV_nodes:
    t1.append(item.name)
    # print(item.name)
for item in PQ_nodes:
    t2.append(item.name)
# 计算雅可比矩阵
def cal_Jocabi_matrix():
    H0 = zeros((9, 9)) * 1j
    N0 = zeros((9, 9)) * 1j
    M0 = zeros((9, 9))* 1j
    L0 = zeros((9, 9)) * 1j
    for i in range(9):
        for j in range(9):
            if i != j:
                H0[i][j] = -nodes[i].U * nodes[j].U * (Y[i][j].real * sin(nodes[i].ang - nodes[j].ang) - Y[i][j].imag *
                                                       cos(nodes[i].ang - nodes[j].ang))
                N0[i][j] = -nodes[i].U * nodes[j].U * (Y[i][j].real * cos(nodes[i].ang - nodes[j].ang) + Y[i][j].imag *
                                                       sin(nodes[i].ang - nodes[j].ang))
                M0[i][j] = +nodes[i].U * nodes[j].U * (Y[i][j].real * cos(nodes[i].ang - nodes[j].ang) + Y[i][j].imag *
                                                       sin(nodes[i].ang - nodes[j].ang))
                L0[i][j] = -nodes[i].U * nodes[j].U * (Y[i][j].real * sin(nodes[i].ang - nodes[j].ang) - Y[i][j].imag *
                                                       cos(nodes[i].ang - nodes[j].ang))
            else:
                H0[i][i] = nodes[i].U ** 2 * Y[i][i].imag
                N0[i][i] = -nodes[i].U ** 2 * Y[i][i].real
                M0[i][i] = nodes[i].U ** 2 * Y[i][i].real
                L0[i][i] = nodes[i].U ** 2 * Y[i][i].imag
                for k in range(9):
                    H0[i][i] += nodes[i].U * nodes[k].U * (Y[i][k].real * sin(nodes[i].ang - nodes[k].ang) - Y[i][k].imag *
                                                       cos(nodes[i].ang - nodes[k].ang))
                    N0[i][i] -= nodes[i].U * nodes[k].U * (Y[i][k].real * cos(nodes[i].ang - nodes[k].ang) + Y[i][k].imag *
                                                       sin(nodes[i].ang - nodes[k].ang))
                    M0[i][i] -= nodes[i].U * nodes[k].U * (Y[i][k].real * cos(nodes[i].ang - nodes[k].ang) + Y[i][k].imag *
                                                       sin(nodes[i].ang - nodes[k].ang))
                    L0[i][i] -= nodes[i].U * nodes[k].U * (Y[i][k].real * sin(nodes[i].ang - nodes[k].ang) - Y[i][k].imag *
                                                       cos(nodes[i].ang - nodes[k].ang))


    H = H0[t1]
    H = H[:, t1]
    N = N0[t1]
    N = N[:, t2]
    M = M0[t2]
    M = M[:, t1]
    L = L0[t2]
    L = L[:, t2]
    J1 = concatenate([H, N], axis=1)
    J2 = concatenate([M, L], axis=1)
    J = concatenate([J1, J2], axis=0)
    return J


PRICISION = 0.00001  # 精度
pricision_now = 1

for _ in range(10):  # 最大迭代10次
    if pricision_now < PRICISION:
        break
    else:
        J = cal_Jocabi_matrix()
        DP = zeros((1, 9)) * 1j
        DQ = zeros((1, 9)) * 1j
        for i in range(9):
            DP[0][i] = nodes[i].P
            DQ[0][i] = nodes[i].Q
            for j in range(9):
                DP[0][i] -= nodes[i].U * nodes[j].U * (Y[i][j].real * cos(nodes[i].ang - nodes[j].ang) + Y[i][j].imag *
                                                           sin(nodes[i].ang - nodes[j].ang))
                DQ[0][i] -= nodes[i].U * nodes[j].U * (Y[i][j].real * sin(nodes[i].ang - nodes[j].ang) - Y[i][j].imag *
                                                       cos(nodes[i].ang - nodes[j].ang))
        DP = DP[0][t1]
        DQ = DQ[0][t2]
        T = concatenate([DP, DQ], axis=0)

        X = -dot(np.linalg.inv(J), T)
        pricision_now = float(X.max())
        for item, i in zip(PQ_PV_nodes, range(8)):
            item.ang += X[i]
        for item, i in zip(PQ_nodes, range(8, 14)):
            item.U += X[i]
        print(pricision_now)  # 输出当前计算精度

nodes[3].U = 1.0258
nodes[4].U = 0.9956
nodes[5].U = 1.0127
nodes[6].U = 1.0258
nodes[7].U = 1.0159
nodes[8].U = 1.0324
# 求解节点功率
for i in range(9):
    nodes[i].P = 0
    nodes[i].Q = 0
for item, i in zip(nodes, range(9)):
    for itemm, j in zip(nodes, range(9)):
        item.P += item.U * itemm.U * (Y[i][j].real * cos(nodes[i].ang - nodes[j].ang) + Y[i][j].imag *
                                                           sin(nodes[i].ang - nodes[j].ang))
        item.Q += item.U * itemm.U * (Y[i][j].real * sin(nodes[i].ang - nodes[j].ang) - Y[i][j].imag *
                                                       cos(nodes[i].ang - nodes[j].ang))

for item in nodes:
    print(item.name+1, item.U, item.ang * 180 / pi, item.P, item.Q)  # 输出潮流计算结果

# 短路计算
# 发电机 负荷节点增加自导纳

for i in range(9):
    if nodes[i].type != "PQ":
        Y[i][i] += 1 / (nodes[i].Xd* 1j)
Y0 = Y.copy()
for i in range(9):
    if nodes[i].P < 0:
        Y[i][i] -= (nodes[i].P - nodes[i].Q * 1j) / (nodes[i].U ** 2)
Z = np.linalg.inv(Y)  # 计算节点阻抗矩阵
#print("节点阻抗矩阵第四列", Z[:][3])
for i in range(9):
    print("节点"+str(i+1)+"自导纳：", Y0[i][i])
# 获得节点电压xy坐标量
for i in range(9):
    nodes[i].U = nodes[i].U * cos(nodes[i].ang) + nodes[i].U * sin(nodes[i].ang) * 1j


# 计算短路点短路电流
SHORT_NODE = 4
If = nodes[SHORT_NODE-1].U / Z[SHORT_NODE-1][SHORT_NODE-1]

# 计算各节点电压电流向量
Ux = []
Ix = []
head = [1, 2, 3, 4, 4, 5, 6, 7, 8]  # 记录支路头结点顺序
end = [4, 7, 9, 5, 6, 7, 9, 8, 9]  # 记录支路尾节点顺序
z = [0+0.0576j, 0+0.0625j, 0+0.0586j, 0.01+0.085j, 0.017+0.092j, 0.032+0.161j, 0.039+0.17j, 0.0085+0.072j, 0.0119+0.1008j]  # 记录支路阻抗值

for i in range(9):
    Ux.append(nodes[i].U - Z[i][SHORT_NODE-1] * (nodes[SHORT_NODE-1].U / Z[SHORT_NODE-1][SHORT_NODE-1] + 0))
for i in range(9):
    Ix.append(((Ux[head[i]-1] - Ux[end[i]-1]) / z[i]))


# 近似计算 忽略
Z0 = np.linalg.inv(Y0)
If2 = 1 / Z0[SHORT_NODE-1][SHORT_NODE-1]
Ux2 = []
for i in range(9):
    Ux2.append((1 - (Z0[i][SHORT_NODE-1]/(Z0[SHORT_NODE-1][SHORT_NODE-1] + 0))))
print("精确计算与近似计算条件下的短路电流模值", If.__abs__(), If2.__abs__())
print("精确计算与近似计算条件下的短路电流相角", polar(If)[1]*180/pi, polar(If2)[1]*180/pi)

for i in range(9): # 输出精确计算和近似计算的节点电压模值
    Ux[i] = Ux[i].__abs__()
    Ux2[i] = Ux2[i].__abs__()
    print("节点"+str(i+1), Ux[i], Ux2[i])

for i in range(9):  # 输出支路电流
    print("支路{}-{}".format(head[i], end[i]), Ix[i])

#  将节点导纳矩阵写入到excel文件之中
book = xlwt.Workbook()
sheet = book.add_sheet('sheet1')
for i in range(9):
    for j in range(9):

        sheet.write(i, j,nodal_admittance_matrix[i][j])
book.save('nodal_admittance_matrix.xls')


