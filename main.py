import convection_diffusion.diff_con as solver

'''

        FINITE VOLUME METHOD for convection_diffusion equation
        точное решение удовлетворяет принципу максимума 0<=U(x)<=1 \forall x \in (0,1)
        при больших Пекле около х=1 имеется пограничный слой

        Dirichlet problem:

        Pe* (Du/Dx) - (D2u/Dx2) = 0 in (0,1)

        u(0) = 0
        u(1) = 1

        --------------------------

        real solution:

        u(x) = (exp(Pe*x)-1)/(exp(Pe) - 1),  du/dx = Pe*(exp(Pe*x))/(exp(Pe)-1)


        Pe = 0, 0.5, 1, 10, 100
'''

if __name__ == "__main__":
    Pe_values = [0.001, 0.5, 1, 10, 100]
    N_values = [11, 21, 41, 81, 161, 641]
    solver.draw_solution(Pe_values, N_values, cols=5)
    solver.draw_error(Pe_values, N_values)

