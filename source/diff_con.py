import numpy as np
from sys import stderr
import time
import thomas
import matplotlib.pyplot as plt


def real_solution_der(x, Pe):
    return Pe * (np.exp(Pe * x)) / (np.exp(Pe) - 1.0)


def real_solution(x, Pe):
    return (np.exp(Pe * x) - 1.0) / (np.exp(Pe) - 1.0)


def thomas_solver(n, d, du, dl, b):
    if d[0] == 0:
        stderr.write('condition w[i]==0 not met')
        exit(-1)
    q, g, u = [_ for _ in np.zeros(n)]
    q[0] = du[0] / d[0]
    g[0] = b[0] / d[0]
    for i in range(1, n):
        w = d[i] - dl[i] * q[i - 1]
        if w == 0:
            stderr.write('condition w[i]==0 not met')
            exit(-1)
        if i != n - 1:
            q[i] = du[i] / w
        g[i] = (b[i] - dl[i] * g[i - 1]) / w

    u[n - 1] = g[n - 1]
    for i in range(n - 2, -1, -1):  # i = N-2, N-1, ... , 0
        u[i] = g[i] - q[i] * u[i + 1]

    return u


def norm_L2(vector1, vector2):
    if len(vector1) != len(vector2):
        print('not equal lengths in norm calculating')
        return 0
    norm = 0
    for i in range(len(vector2)):
        norm += ((vector2[i] - vector1[i]) ** 2)
    return np.sqrt(norm)


def counterflow(N, Pe=10):
    dirichlet0 = 0
    dirichletN = 1
    h = 1 / (N)
    # N = N - 1
    d = np.zeros(N)
    du = np.zeros(N)
    dl = np.zeros(N)
    b = np.zeros(N)
    for i in range(0, N):
        b[i] = 0
        du[i] = - 1 / (h ** 2)
        d[i] = Pe / h + 2 / (h ** 2)
        dl[i] = -Pe / h - 1 / (h ** 2)
    b[0] = - dirichlet0 * (-Pe / h - 1 / (h ** 2))
    b[N - 1] = - dirichletN * (- 1 / (h ** 2))

    # print('d',d)
    # print('du',du)
    # print('dl', dl)
    # print('b',b)
    return d, du, dl, b


# between scheme
def solver_CD(N, Pe):
    h = 1. / (N - 1)

    boundary_left = 0
    boundary_right = 1.
    N = N - 1
    d = np.ones(N) * 2 / (h ** 2)  # diagonal
    dl = np.ones(N) * (-1 * Pe / (2 * h) - 1 / (h ** 2))  # upper diagonal
    du = np.ones(N) * (1 * Pe / (2 * h) - 1 / (h ** 2))  # lower diagonal

    b = np.zeros(N)

    b[0] = b[0] - dl[0] * boundary_left
    b[N - 1] = b[N - 1] - du[N - 1] * boundary_right

    A = thomas.ThreeDiagMatrix(d, du, dl)
    # print('du: ', du)
    # print('d: ', d)
    # print('dl: ', dl)
    # print('b: ', b)
    u_first_method = A.thomas_solver(b)
    real_sol = np.zeros(N)
    x = []
    for i in range(N):
        x.append(i * h + h / 2)
        # print('x: ',x[i])
        real_sol[i] = real_solution(x[i], Pe)

    # print('numeric: ', u_first_method)
    # print('real:    ', real_sol)

    print('----------------------------------------------------------')
    print('Pe', Pe, ' N = ', N + 1, ' L2 norm: ', norm_L2(u_first_method, real_sol))


    plt.plot(x, u_first_method, label='numeric')
    plt.plot(x, real_sol, label='real')
    plt.legend()

    where_to_save = f'../images/Pe_{Pe}_N_{N+1}.png'
    plt.savefig(where_to_save)
    plt.close()

    return u_first_method



def testAli(N, Pe=10):
    h = 1 / (N - 1)
    N = N - 1
    d, du, dl, b = counterflow(N, Pe)
    x = thomas_solver(N, d, du, dl, b)
    analitic = []
    for i in range(0, N):
        analitic.append(real_solution(i * h + h / 2, Pe))
    res = norm_L2(analitic, x)
    print(f'N:{N}  norma: {res}')


def solve_for_pe(Pe):
    N = 11
    solver_CD(N, Pe)

    N = 21
    solver_CD(N, Pe)

    N = 41
    solver_CD(N, Pe)

    N = 81
    solver_CD(N, Pe)

    N = 161
    solver_CD(N, Pe)

    N = 641
    solver_CD(N, Pe)



if __name__ == "__main__":
    Pe_values = [0.001, 0.5, 1, 10, 100]
    N_values = [11, 21, 41, 81, 161, 641]
    for Pe in Pe_values:
        for N in N_values:
            solver_CD(N, Pe)


