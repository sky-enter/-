import numpy as np
import pandas as pd
import numpy.linalg as lg

# data = pd.read_excel('生态大楼K.xls')
data = pd.read_csv('F:\\解超F\\数据.csv')
print(data)

class Dynamic(object):
    def __init__(self, heat_flow, inside, outside, number_equations=30, det_t=0.5, m=3, r=5):
        self.M = number_equations  # 方程的个数
        self.det = det_t  # 时间间隔h
        self.N = len(heat_flow) - 2  # 数据总数
        self.P = self.N - self.M  # 用于求和的P组数据
        self.m = m  # 时间常数的个数不大于3
        self.r = r  # 时间常数的不变比率在3-10之间
        self.tao1 = self.det / 10
        self.tt = []
        self.s = []
        ss = 100000
        rr = 0
        tao = 0
        while self.tao1 <= self.P * self.det / 2:
            # while True:
            summ, R = self.sum_s2(heat_flow, inside, outside)
            if summ < 20000:
                self.s.append(summ)
                self.tt.append(self.tao1)
            if summ < ss:
                ss = summ
                rr = R
                tao = self.tao1
            # print('时间常数=%.2f,S2=%.1f,R=%.3f' % (self.tao1, summ, R))
            self.tao1 = self.det / 10 + self.tao1
        print('时间常数=%.2f,S2=%.1f,R=%.3f' % (tao, ss, rr))

    def matrix_x(self, heat_flow, inside, outside):
        heat = heat_flow
        inst = inside
        outs = outside
        heat_flow = heat_flow[1:-1]
        inside = inside[1:-1]
        outside = outside[1:-1]

        xp = []
        tao1 = self.tao1
        tao = [tao1, tao1 / self.r, tao1 / np.sqrt(self.r)]

        q1 = heat_flow[self.P: self.N]  # 最后M个热流密度数据
        for i in range(self.P, self.N):

            x = []
            for j in range(1, 2 * self.m + 4):
                if j == 1:
                    x.append(inside[i] - outside[i])  # 第i行第1列：室内外温差
                elif j == 2:
                    x.append((inside[i] - inside[i - 1]) / self.det)  # 第i行第2列：内表面温度导数
                elif j == 3:
                    x.append((outside[i] - outside[i - 1]) / self.det)  # 第i行第3列：外表面温度导数
                else:
                    if j % 2 == 0:  # 偶数列的时间行数的选择
                        tt = (j - 2) / 2 - 1
                    else:  # 奇数列的时间行数的选择
                        tt = (j - 3) / 2 - 1
                    bata = np.exp(-self.det / tao[int(tt)])  # 时间常数的指数函数
                    xe = 0
                    for k in range(i - self.P, i):
                        if j % 2 == 0:  # 偶数列的计算

                            xee = (inside[k] - inside[k - 1]) / self.det
                            xe = xe + xee * (1 - bata) * bata * (i - k)
                            if k == 0:
                                xee = (inside[k] - inst[k]) / self.det
                                xe = xe + xee * (1 - bata) * bata * (i - k)
                        else:
                            xee = (outside[k] - outside[k - 1]) / self.det
                            xe = xe + xee * (1 - bata) * bata * (i - k)
                            if k == 0:
                                xee = (outside[k] - outs[k]) / self.det
                                xe = xe + xee * (1 - bata) * bata * (i - k)
                    x.append(xe)
            xp.append(x)
        xp = np.array(xp)
        q1 = np.array(q1)
        return xp, q1

    def sum_s2(self, heat_flow, inside, outside):
        x, q = self.matrix_x(heat_flow, inside, outside)
        Y = x.T.dot(x)
        try:
            Y = lg.inv(Y)
        except:
            return 100000, 0
        # else:
        # print(Y)

        a3 = Y.dot(x.T)
        zz = a3.dot(q.T)
        qz = x.dot(zz)
        s_2 = (q - qz) ** 2
        sum_s_2 = np.sum(s_2)

        self.I = np.sqrt(sum_s_2 * Y[0, 0] / (self.M - 2 * self.m - 4))

        return sum_s_2, 1 / zz[0]

# r = Dynamic(q[0:72], ti[0:72], te[0:72])
# r = Dynamic(q[44:188], ti[44:188], te[44:188])
# print(r.I)
# import matplotlib.pyplot as plt
#
# plt.plot(r.tt, r.s)
# # plt.plot(r.tt, r.R)
#
# plt.show()
