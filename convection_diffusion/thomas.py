import random
import time

import numpy as np
import scipy as sp
from sys import stderr

'''---------------------------------
    forward
    q1 = b1/a1
    g1 = d1/a1
    for 2,..,N
        wi = ai -ci*q[i-1] !=0
        qi = bi/wi 
        gi = (di-ci*g[i-1])/wi
    U_n = g_n
    backward
    for i = n-1, ... , 1
        ui=gi-qi*u[i+1]
-----------------------------------'''


class ThreeDiagMatrix:
    def __init__(self, d, du, dl):
        self.d = d
        self.du = du
        self.dl = dl

    # def __init__(self, N):
    #     self.d = np.zeros(N)
    #     self.du = np.zeros(N)
    #     self.dl = np.zeros(N)

    def fill_random(self, N):
        for i in range(0, n):
            self.d[i] = int(random.random() * 100)
            if i > 0:
                self.dl[i] = int(random.random() * 100)
            if i < n - 1:
                self.du[i] = int(random.random() * 100)

    def thomas_solver(self, d):
        if self.d[0] == 0:
            print(i)
            stderr.write('d[0]==0')
            exit(-1)
        N = len(self.d)

        q = np.zeros(N)
        g = np.zeros(N)
        # w = np.zeros(N)
        u = np.zeros(N)
        start_time = time.time()
        # forward
        q[0] = self.du[0] / self.d[0]
        g[0] = d[0] / self.d[0]
        for i in range(1, N):
            w = self.d[i] - self.dl[i] * q[i - 1]
            # print('i: ',i,' w: ',w)
            # print('d[i] = ', self.d[i], '  dl[i] = ', self.dl[i], ' q[i-1] = ', q[i - 1])
            if w == 0:
                print('index: ', i)
                print('d[i] = ', self.d[i], '  dl[i] = ', self.dl[i], ' q[i-1] = ', q[i - 1])
                print('condition w[i]==0 not met')
                exit(-1)
            if i != N - 1:
                q[i] = self.du[i] / w
            g[i] = (d[i] - self.dl[i] * g[i - 1]) / w

        # backward
        u[N - 1] = g[N - 1]
        for i in range(N - 2, -1, -1):  # i = N-2, N-1, ... , 0
            u[i] = g[i] - q[i] * u[i + 1]
        end_time = time.time()
        # print('time of my solver: ', (time.time() - start_time) / 1000, 'ms')
        # print('L2 norm: ', sp.linalg.norm(self.mul_mat_vec(u) - d, 2))
        return u

    def mul_mat_vec(self, b):
        if len(self.d) != len(b):
            exit(-1)
        res = np.zeros(len(self.d))
        for i in range(len(self.d)):
            if i > 0:
                res[i] += self.dl[i] * b[i - 1]
            res[i] += self.d[i] * b[i]
            if i < len(self.d) - 1:
                res[i] += self.du[i] * b[i + 1]
        return res

#
# if __name__ == '__main__':
#     N = 100
#     d = np.ones(N) * 2
#     du = np.ones(N) * (-1)
#     dl = np.ones(N) * (-1)
#     b = np.ones(N)
#
#     A = ThreeDiagMatrix(d, du, dl)
#     print('------------------------')
#     # print('du: ', A.du)
#     # print('d: ', A.d)
#     # print('dl: ', A.dl)
#
#     norm = 0
#     my_solution = A.thomas_solver(b)
#     res = A.mul_mat_vec(my_solution)
#     for i in range(len(b)):
#         norm += ((res[i] - b[i]) * (res[i] - b[i]))
#
#     print('L2-norm: ', np.sqrt(norm))
#
#     N = 100
#     d = 2 + np.random.rand(N)
#     d = -np.random.rand(N)
#     du = -np.random.rand(N)
#     dl = -np.random.rand(N)
#     b = 2*np.random.rand(N)-1
#
#     A = ThreeDiagMatrix(d, du, dl)
#     print('------------------------')
#     # print('du: ', A.du)
#     # print('d: ', A.d)
#     # print('dl: ', A.dl)
#
#     norm = 0
#     my_solution = A.thomas_solver(b)
#     res = A.mul_mat_vec(my_solution)
#     for i in range(len(b)):
#         norm += ((res[i] - b[i]) * (res[i] - b[i]))
#
#     print('L2-norm: ', np.sqrt(norm))
