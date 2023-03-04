import numpy as np
from sys import stderr
import time
import convection_diffusion.thomas as thomas
import matplotlib.pyplot as plt
from typing import Tuple, List


def integral(left, right, Pe, h):
    return (((np.exp(Pe * right) / Pe - right) / (np.exp(Pe) - 1)) - (
            (np.exp(Pe * left) / Pe - left) / (np.exp(Pe) - 1))) / h


def real_solution_der(x: float, Pe: float):
    return Pe * (np.exp(Pe * x)) / (np.exp(Pe) - 1.0)


def real_solution(x: float, Pe: float):
    return (np.exp(Pe * x) - 1.0) / (np.exp(Pe) - 1.0)


def norm_L2(vector1: np.array, vector2: np.array, h) -> float:
    if len(vector1) != len(vector2):
        print('not equal lengths in norm calculating')
        return 0
    norm = 0
    for i in range(len(vector2)):
        norm += ((vector2[i] - vector1[i]) ** 2)
    return np.sqrt(norm * h)


def counterflow(N: int, Pe: float = 10) -> Tuple[list, np.array, np.array]:
    dirichlet0, dirichletN = 0, 1
    h = 1 / (N - 1)
    N = N - 1
    d = np.zeros(N)
    du = np.zeros(N)
    dl = np.zeros(N)
    b = np.zeros(N)

    b[0] = - dirichlet0 * (-2 * Pe / h - 2 / (h ** 2))

    du[0] = - 1 / (h ** 2)
    d[0] = 2 * Pe / h + 3 / (h ** 2)
    dl[0] = - 2 * Pe / h - 2 / (h ** 2)
    for i in range(1, N - 1):
        b[i] = 0
        du[i] = - 1 / (h ** 2)
        d[i] = Pe / h + 2 / (h ** 2)
        dl[i] = -Pe / h - 1 / (h ** 2)

    b[N - 1] = - dirichletN * (- 2 / (h ** 2))

    du[N - 1] = - 2 / (h ** 2)
    d[N - 1] = Pe / h + 3 / (h ** 2)
    dl[N - 1] = - Pe / h - 1 / (h ** 2)

    A = thomas.ThreeDiagMatrix(d, du, dl)
    u = A.thomas_solver(b)
    analitic = np.zeros(N)
    x = np.zeros(N)
    for i in range(N):
        analitic[i] = integral(i * h, (i + 1) * h, Pe, h)
        x[i] = (i * h + h / 2)

    print(f'Pe {Pe} N = {N + 1} L2 norm: {norm_L2(u, analitic, h)}')

    return x, u, analitic


# between scheme
def solver_CD(N: int, Pe: float) -> Tuple[list, np.array, np.array]:
    h = 1. / (N - 1)

    boundary_left = 0
    boundary_right = 1.
    N = N - 1
    d = np.ones(N) * 2 / (h ** 2)  # diagonal
    dl = np.ones(N) * (-1 * Pe / (2 * h) - 1 / (h ** 2))  # upper diagonal
    du = np.ones(N) * (1 * Pe / (2 * h) - 1 / (h ** 2))  # lower diagonal

    b = np.zeros(N)
    dl[0] = -Pe / h - 2 / (h ** 2)
    d[0] = Pe / (2 * h) + 3 / (h ** 2)
    du[0] = Pe / (2 * h) - 1 / (h ** 2)

    dl[N - 1] = - Pe / (2 * h) - 1 / (h ** 2)
    d[N - 1] = Pe / (2 * h) + 3 / (h ** 2)
    du[N - 1] = Pe / h - 2 / (h ** 2)
    b[0] = b[0] - dl[0] * boundary_left
    b[N - 1] = b[N - 1] - du[N - 1] * boundary_right

    A = thomas.ThreeDiagMatrix(d, du, dl)
    u_first_method = A.thomas_solver(b)
    real_sol = np.zeros(N)
    x = np.zeros(N)
    for i in range(N):
        real_sol[i] = integral(i * h, (i + 1) * h, Pe, h)
        x[i] = (i * h + h / 2)

    print(f'Pe {Pe} N = {N + 1} L2 norm: {norm_L2(u_first_method, real_sol, h)}')

    return x, u_first_method, real_sol


def draw_error(Pe_values: List[float], N_values: List[int]) -> None:
    rows, cols = 2, len(Pe_values)
    fig, axs = plt.subplots(rows, cols, figsize=(20, 7), constrained_layout=True)

    for i, Pe in enumerate(Pe_values):

        array_l2_CD = []
        array_l2_BD = []
        for N in N_values:
            h = 1. / (N - 1)

            xcd, ucd, solcd = solver_CD(N, Pe)  # solver CD
            array_l2_CD.append(norm_L2(ucd, solcd, h))
            xbd, ubd, solbd = counterflow(N, Pe)  # solver BD
            array_l2_BD.append(norm_L2(ubd, solbd, h))

        axs[0, i].plot(np.log(N_values), np.log(array_l2_CD), label=f"Pe = {Pe}", c='green')
        axs[0, i].grid(True)
        axs[0, i].legend()

        axs[1, i].plot(np.log(N_values), np.log(array_l2_BD), label=f"Pe = {Pe}", c='red')
        axs[1, i].grid(True)
        axs[1, i].legend()

    axs[0, 2].set_title(" Norm L2 Main value between cells")

    axs[1, 2].set_title(" Norm L2 Anti-Flow Scheme")

    path_to_save = f'images/convection-diffusion/L2_norm.png'
    plt.savefig(path_to_save)
    plt.close()


def draw_solution(Pe_values: List[float], N_values: List[int], cols: int = 5) -> None:
    """
    :param Pe_values: list of Pe values;
    :param N_values: list of N values;
    :param cols: 2 or 5 (default);
    :return: figures;
    """

    if cols == 2:
        ROWS = (len(N_values) // 2 + len(N_values) % 2) * 2
    else:
        ROWS, cols = (len(Pe_values) // 5 + len(Pe_values) % 5) * 2, 5

    axsIndexes = [(i, j) for i in range(ROWS) for j in range(cols)]

    for N in N_values:
        print("------------------------------------------")
        if cols == 2:
            fig, axs = plt.subplots(ROWS, cols, figsize=(10, 15), constrained_layout=True)
        else:
            fig, axs = plt.subplots(ROWS, cols, figsize=(20, 7), constrained_layout=True)
        plt.suptitle(f"Steps = {N - 1}", fontsize=20)
        for i, Pe in enumerate(Pe_values):
            x, u, sol = solver_CD(N, Pe)  # solver CD
            axs[axsIndexes[i]].plot(x, sol, label=f"Real Solution, Pe = {Pe}", c='green')
            axs[axsIndexes[i]].plot(x, u, label=f"Numeric Solution, Pe = {Pe}", c='red')

            x, u, sol = counterflow(N, Pe)  # solver BD
            axs[axsIndexes[i + len(Pe_values)]].plot(x, sol, label=f"Real Solution, Pe = {Pe}", c='green')
            axs[axsIndexes[i + len(Pe_values)]].plot(x, u, label=f"Numeric Solution, Pe = {Pe}", c='red')

            axs[axsIndexes[i]].set_title("Main value between cells")
            axs[axsIndexes[i]].grid(True)
            axs[axsIndexes[i]].legend()

            axs[axsIndexes[i + len(Pe_values)]].set_title("Anti-Flow Scheme")
            axs[axsIndexes[i + len(Pe_values)]].grid(True)
            axs[axsIndexes[i + len(Pe_values)]].legend()

        if not (len(N_values) % cols):
            fig.delaxes(axs[-1][-1])

        path_to_save = f'images/convection-diffusion/N_{N - 1}.png'
        plt.savefig(path_to_save)
        plt.close()
