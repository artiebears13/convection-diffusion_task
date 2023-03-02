import numpy as np
from sys import stderr
import time


def real_solution_der(x, Pe = 1):
    return Pe * (np.exp(Pe * x )) / (np.exp(Pe) - 1.0)


def real_solution(x, Pe = 1):
    return (np.exp(Pe * x) - 1.0) / (np.exp(Pe) - 1.0)

def thomas_solver(n, d, du, dl, b):
    if d[0] == 0:
        stderr.write('condition w[i]==0 not met')
        exit(-1)
    q = np.zeros(n)
    g = np.zeros(n)

    u = np.zeros(n)
    q[0] = du[0]/d[0]
    g[0] = b[0] / d[0]
    for i in range(1, n):
        w = d[i] - dl[i] * q[i - 1]
        if w == 0:
            stderr.write('condition w[i]==0 not met')
            exit(-1)
        if i != n-1:
            q[i] = du[i]/w
        g[i] = (b[i] - dl[i] * g[i-1]) / w

    u[n - 1] = g[n - 1]
    for i in range(n - 2, -1, -1):  # i = N-2, N-1, ... , 0
        u[i] = g[i] - q[i] * u[i + 1]

    return u



def main():




if __name__ == '__main__':
    main()
