from plot_drawer import show_3d_graph
from stat_solver import *
from data_former import *
from particle import Particle
import numpy as np


def f(x, y, x0, y0):
    alpha = 0.1
    res = 10 * np.exp(-alpha * ((x - x0)**2 + (y-y0)**2))
    return res

def f_by_index(i, j, X, Y, x0, y0):
    return f(X[i], Y[j], x0, y0)


def main():
    a, b, c, d = get_data_from_file()
    X, Y, Z = form_plane(a, b)

    lenX = len(X)
    lenY = len(Y)

    u0 = 300
    form_circuit(Z, u0)

    in_area = form_in_area(a, b, c, d, lenX, lenY)
    define_in_area(Z, in_area, X[1] - X[0], Y[1] - Y[0])
    #define_output_area(Z, in_area, X[1] - X[0], Y[1] - Y[0], u0)

    
    # Точка нагрева
    x0 = (3 / 2) * in_area['border'][0] + in_area['size'][0]
    y0 = (1 / 2) * in_area['border'][1]
    heat_point = [x0, y0]
    print('(X0, Y0): (', x0, ';', y0, ')')


    # Стохастический метод
    M = 10
    for i in range(lenX):
        for j in range(lenY):
            stochastic_solve_point(i, j, X, Y, Z, heat_point, in_area, u0, CheckBound, rightPart, M)

    show_3d_graph(X, Y, Z, center=(heat_point[0], heat_point[1]))


if __name__ == "__main__":
    main()