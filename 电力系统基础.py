import pandas as pd
from math import sqrt, atan
import xlwt
import datetime

"""
----------------------------------------------------------------------------------------------------------------------------
作者：李予辉 电气二班 学号： 3017234232
课程名称：电力系统基础
完成时间：2019.12.16
----------------------------------------------------------------------------------------------------------------------------
"""
def generate_result_file(status, node_num, num):  # status为检测输入变量是否符合规则 num 迭代次数
    """用于生成结果文件的函数"""
    workbook = xlwt.Workbook(encoding='utf-8')
    worksheet = workbook.add_sheet('NodesResult')
    for i in range(20):
        worksheet.col(i).width = 5000
    if times > 500:
        status = '迭代次数超过500次，可能不收敛，停止迭代'
    # 写入数据
    # 写入输出状态信息
    temp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')  # 现在时间

    worksheet.write_merge(0, 0, 0, 10, "结果文件输出成功，输出时间：{}. 输出状态：{}".format(temp, status))
    worksheet.write(1, 0, '节点')
    worksheet.write(1, 1, '父节点')
    worksheet.write(1, 2, '子节点')
    worksheet.write(1, 3, '节点电压kV')
    worksheet.write(1, 4, '电压相角deg')
    worksheet.write(1, 5, '前向支路始端有功MW')
    worksheet.write(1, 6, '前向支路始端无功MVar')
    worksheet.write(1, 7, '前向支路末端有功MW')
    worksheet.write(1, 8, '前向支路末端无功MVar')
    worksheet.write(1, 9, '前向支路电压损耗kV')
    for i in range(node_num):
        for j in range(10):
            worksheet.write(i+2, j, output1[i][j])

    worksheet.write(node_num + 2, 0, '全网的总电源有功MW')
    worksheet.write(node_num + 2, 1, '总负荷有功MW')
    worksheet.write(node_num + 2, 2, '全网有功损耗MW')
    worksheet.write(node_num + 2, 3, '网损率')
    worksheet.write(node_num + 2, 4, '最低电压')
    worksheet.write(node_num + 2, 5, '最低电压点')
    worksheet.write(node_num + 2, 6, '迭代总次数')
    for i in range(7):
        worksheet.write(node_num + 3, i, output2[i])

    for i in range(num):
        worksheet.write((i+2)*(node_num+2)-1, 0, '迭代次数{}'.format(str(i+1)))
        worksheet.write((i+2)*(node_num+2), 0, '节点')
        worksheet.write((i+2)*(node_num+2), 1, '节点电压kV')
        worksheet.write((i+2)*(node_num+2), 2, '节点负荷有功MW')
        worksheet.write((i+2)*(node_num+2), 3, '节点负荷无功MVar')
        for j in range(node_num):
            for k in range(4):
                worksheet.write((i + 2) * (node_num + 2)+1+j, k, output_times[i][j][k])


    # 写入其他参数

    workbook.save("数据输入输出/ResultOutput.xls")


class Node():
    def __init__(self, P_node, C_node, L, rl, xl, b0, Pload, Qload, k, Sn, P0, Pk, Uk010, I0010, RT, XT, PT_0, QT_0, UN, U0):
        self.P_node = P_node  # 父节点
        self.C_node = C_node  # 子节点
        self.L = L  # 线路长度
        self.rl = rl
        self.xl = xl
        self.b0 = b0
        self.Pload = Pload
        self.Qload = Qload
        self.k = k
        self.cal_k()
        self.Sn = Sn
        self.P0 = P0
        self.Pk = Pk
        self.Uk010 = Uk010
        self.I0010 = I0010
        self.RT_temp = RT
        self.XT_temp = XT
        self.PT_0 = PT_0
        self.QT_0 = QT_0
        self.UN = UN  # 线路额定电压
        self.U = self.UN * self.k[1] / self.k[0]

        self.RT = self.cal_TRT()
        self.XT = self.cal_TXT()
        self.S0_P = self.cal_TS0_P()  # 计算变压器空载损耗
        self.S0_Q = self.cal_TS0_Q()


        self.S1_P = 0
        self.S2_Q = 0
        self.S1_Q = 0
        self.S2_P = 0
        self.S_P = 0
        self.S_Q = 0
        self.child_P = 0
        self.child_Q = 0
        self.U0 = U0
        self.U1 = 0
        self.U2 = 0
        self.big_delta_U = 0
        self.small_delta_U = 0
        #print(self.RT, self.XT)

    def cal_k(self):
        # 计算变压器变比
        temp = self.k.split('，')
        self.k = []
        self.k.append(int(temp[0]))
        self.k.append(int(temp[1]))

    def change_P_C_nodes(self):
        if self.P_node != "Null":
            self.P_node = node_object[int(self.P_node) - 1]
        else:
            self.P_node = None
            self.U = self.U0
        if self.C_node != "Null":
            temp = str(self.C_node)
            temp = temp.split('，')
            self.C_node = []
            for a in temp:
                self.C_node.append(node_object[int(a) - 1])
        else:
            self.C_node = []

    def cal_TRT(self):
        if self.RT_temp != 0:
            return self.RT_temp
        else:
            if self.Sn == 0:
                return 0
            else:
                return self.Pk * self.UN * self.UN / self.Sn/self.Sn * 1000

    def cal_TXT(self):
        if self.XT_temp != 0:
            return self.XT_temp
        else:
            if self.Sn == 0:
                return 0
            else:
                return self.Uk010 * self.UN * self.UN / 100 / self.Sn * 1000

    def cal_TS0_P(self):
         if self.PT_0 != 0:
             return self.PT_0
         else:
             if self.Sn == 0:
                 return 0
             else:
                return self.P0

    def cal_TS0_Q(self):
        if self.QT_0 != 0:
            return self.QT_0
        else:
            if self.Sn == 0:
                return 0
            else:
                return self.I0010 * self.Sn / 100

    def cal_node_S(self):  # 从后向前推功率
        self.child_P = 0
        self.child_Q = 0
        if len(self.C_node) != 0:
            for child in self.C_node:
                child.cal_node_S()
                self.child_P += child.S_P
                self.child_Q += child.S_Q
        delta_QB = -0.5 * self.UN * self.UN * self.L * self.b0
        self.S2_P = self.Pload + self.child_P
        self.S2_Q = self.Qload + self.child_Q + delta_QB
        self.delta_S_P = (self.S2_P * self.S2_P + self.S2_Q * self.S2_Q) / self.U / (self.k[0]/self.k[1]) / self.U / (self.k[0]/self.k[1]) * (self.L * self.rl + self.RT)
        self.delta_S_Q = (self.S2_P * self.S2_P + self.S2_Q * self.S2_Q) / self.U / (self.k[0]/self.k[1]) / self.U / (self.k[0]/self.k[1]) * (self.L * self.xl + self.XT)
        self.S1_P = self.S2_P + self.delta_S_P + self.S0_P
        self.S1_Q = self.S2_Q + self.delta_S_Q + self.S0_Q
        self.S_P = self.S1_P
        self.S_Q = self.S1_Q + delta_QB
        #print(self.U)
        #print(delta_QB, self.S2_P, self.S2_Q, self.delta_S_P, self.delta_S_Q, self.S1_P, self.S1_Q, self.S_P, self.S_Q)

    def cal_node_U(self):

        if self.P_node is not None:
            self.big_delta_U = (self.S1_P * (self.L*self.rl + self.RT) + self.S1_Q * (self.L*self.xl + self.XT)) / self.P_node.U
            self.small_delta_U = (self.S1_P * (self.L*self.xl + self.XT) - self.S1_Q * (self.L*self.rl + self.RT)) / self.P_node.U
            self.U = (sqrt((self.P_node.U - self.big_delta_U) * (self.P_node.U - self.big_delta_U) + self.small_delta_U * self.small_delta_U)) * (self.k[1]/self.k[0])
            self.big_delta_U2 = (self.S2_P * (self.L*self.rl + self.RT) + self.S2_Q * (self.L*self.xl + self.XT)) / self.U
            self.small_delta_U2 = (self.S2_P * (self.L*self.xl + self.XT) - self.S2_Q * (self.L*self.rl + self.RT)) / self.U
            #print(self.U)
            self.U1 = sqrt((self.U+self.big_delta_U2)**2 + self.small_delta_U2**2)
            self.U2 = sqrt((self.U1 - self.big_delta_U)**2 + self.small_delta_U**2)

        if len(self.C_node) != 0:
            for child in self.C_node:
                child.cal_node_U()


if __name__ == '__main__':

    output_times = []  # 每次迭代的数据
    output_now = []  # 本次迭代的数据
    output_first = []  # 初始数据
    num = 0  # 记录迭代次数
    data = pd.read_excel('数据输入输出/NodeData.xlsx')
    node_num = 0

    for i in range(len(data['节点'].values)):  # 计算节点个数,默认节点1为头结点
        try:
            a = int(data['节点'].values[i])
            if a != 0:
                node_num += 1
        except:
            continue

    node_data = []
    node_S_result = []
    node_U_result = []
    U0 = data['配电网首端电压（kV）'].values[0]
    UN = data['配电网首端电压（kV）'].values[2]
    error = data['配电网首端电压（kV）'].values[4]

    for i in range(node_num):
        temp = (data.values[i][1:][1], data.values[i][1:][2], data.values[i][1:][3], data.values[i][1:][4], data.values[i][1:][5], data.values[i][1:][6], data.values[i][1:][7], data.values[i][1:][8], data.values[i][1:][9], data.values[i][1:][10], data.values[i][1:][11], data.values[i][1:][12], data.values[i][1:][13], data.values[i][1:][14], data.values[i][1:][15], data.values[i][1:][16], data.values[i][1:][17], data.values[i][1:][18])
        node_data.append(temp)

    node_object = []
    for node in node_data:
        temp1 = None
        temp2 = None
        if node[0] is not 'Null':
            temp1 = node[0]
        if node[1] is not 'Null':
            temp2 = node[1]
        node_object.append(Node(temp1, temp2, node[2], node[3], node[4], node[5], node[6], node[7], node[8], node[9] ,node[10] ,node[11], node[12], node[13], node[14], node[15], node[16], node[17], UN, U0))

    for i in node_object:
        i.change_P_C_nodes()
        output_first.append(i.U)

    err = 100
    times = 0
    while err >= error and times <= 500:
        times += 1
        node_object[0].cal_node_S()
        node_object[0].cal_node_U()
        output_now.clear()
        for i in node_object:
            output_now.append(i.U)
        err = 0
        temp2 = []
        for i in range(node_num):
            if abs(output_now[i] - output_first[i]) > err:
                #print(output_now[i], output_first[i])
                err = abs(output_now[i] - output_first[i])
            output_first[i] = output_now[i]
            temp1 = []
            temp1.append(str(i+1))
            temp1.append(node_object[i].U)
            temp1.append(node_object[i].S_P)
            temp1.append(node_object[i].S_Q)
            temp2.append(temp1)

        output_times.append(temp2)
        num += 1
    #print(output_times)
    output1 = []
    output2 = []
    all_power_P = node_object[0].S_P
    all_P = 0
    all_lost = 0
    lowest_U = 10000000
    lowest_node = 0
    for i in range(len(node_object)):
        temp = []
        temp.append(str(i+1))
        temp.append(data.values[i][1:][1])
        temp.append(data.values[i][1:][2])
        temp.append(node_object[i].U)
        temp.append(atan(node_object[i].small_delta_U/(node_object[i].U - node_object[i].big_delta_U))*180/3.14159265358979)
        temp.append(node_object[i].S1_P)
        temp.append(node_object[i].S1_Q)
        temp.append(node_object[i].S2_P)
        temp.append(node_object[i].S2_Q)

        if i==0:
            temp.append(0)
        else:
            temp.append(node_object[i].U1 - node_object[i].U2)

        output1.append(temp)


        all_P += node_object[i].S_P
        all_lost += node_object[i].delta_S_P + node_object[i].S0_P
        if node_object[i].U < lowest_U:
            lowest_U = node_object[i].U
            lowest_node = i+1
    lost_rate = all_lost / all_power_P
    output2.append(all_power_P)
    output2.append(all_P)
    output2.append(all_lost)
    output2.append(lost_rate)
    output2.append(lowest_U)
    output2.append(str(lowest_node))
    output2.append(num)

    for j in range(4):
        output1[0][j+5] = 0

    generate_result_file('正常', node_num, num)
