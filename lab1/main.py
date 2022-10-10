from plot_drawer import show_3d_graph
from stat_solver import *
from data_former import *
from particle import Particle
import numpy as np


k = 0.0134 # 0.001 # 1 #

def f(x, y, x0, y0):
    alpha = 0.1
    return 10 * np.exp(-alpha * ((x - x0)**2 + (y-y0)**2))

def f_by_index(i, j, X, Y, x0, y0):
    return f(X[i], Y[j], x0, y0)

def rightPart(x, y, x0, y0):
    return f(x, y, x0, y0) / k


# краевое 2ого рода
# -k du/dx = F0

# краевое 3его рода
# -k du/dx = a1 (u - u0)

F0 = 30
a0 = 0.1
def CheckBound(p: Particle, X, Y, x0, y0, in_area, u0):
    lenX = len(X)
    lenY = len(Y)

    hX = X[1] - X[0]
    hY = Y[1] - Y[0]

    cX = in_area['border'][0]
    cY = in_area['border'][1]

    indC = in_area['size'][0]
    indD = in_area['size'][1]

    icX = int(cX / hX)
    jcY = int(cY / hY)

    a0 = 0.1

    # левая граница внешнего прямоугольника
    if p.i <= 0:
        # return p.partialReflection(u0, f_by_index, a0, k, hx, di = +1)
        return p.reflection(F0, k, f_by_index, X, Y, x0, y0, hX, di = +1) # краевое 2ого рода
        # return p.absorption(u0) # краевое 1ого рода

    # правая граница внешнего прямоугольника
    if p.i >= len(X) - 1:
        # return p.partialReflection(u0, f_by_index, a0, k, hx, di = -1)
        return p.reflection(F0, k, f_by_index, X, Y, x0, y0, hX, di = -1)
        # return p.absorption(u0)

    # нижняя граница внешнего прямоугольника
    if p.j <= 0:
        # return p.partialReflection(u0, f_by_index, a0, k, hy, dj = +1)
        return p.reflection(-F0, k, f_by_index, X, Y, x0, y0, hY, dj = +1)
        # return p.absorption(u0)

    # верхняя граница внешнего прямоугольника
    if p.j >= len(Y) - 1:
        # return p.partialReflection(u0, f_by_index, a0, k, hy, dj = -1)
        return p.reflection(-F0, k, f_by_index, X, Y, x0, y0, hY, dj = -1)
        # return p.absorption(u0)

    # во внутреннем прямоугольнике
    if (jcY <= p.j and p.j <= jcY + indD) and (icX <= p.i and p.i <= icX + indC):
        # левая граница внутреннего прямоугольника
        if (icX == p.i):
            return p.partialReflection(u0, f_by_index, X, Y, x0, x0, a0, k, hX, di = -1)
            # return p.absorption(u0)

        # правая граница внутреннего прямоугольника
        if (p.i == icX + indC):
            return p.partialReflection(u0, f_by_index, X, Y, x0, y0, a0, k, hX, di = +1)
            # return p.absorption(u0)

        # нижняя граница внутреннего прямоугольника
        if (jcY == p.j):
            return p.partialReflection(u0, f_by_index, X, Y, x0, y0, a0, k, hY, dj = -1)
            # return p.absorption(u0 - 1)

        # верхняя граница внутреннего прямоугольника
        if (p.j == jcY + indD):
            return p.partialReflection(u0, f_by_index, X, Y, x0, y0, a0, k, hY, dj = +1)
            # return p.reflection(F0, k, f_by_index, hy, dj = +1)
            # return p.absorption(u0)

        # отверстие внутреннего прямоугольника
        p.absorption(np.nan)


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
    print('(X0, Y0): (', x0, ';', y0, ')')


    # Стохастический метод
    M = 10
    for i in range(lenX):
        for j in range(lenY):
            stochastic_solve_point(i, j, X, Y, Z, x0, y0, in_area, u0, CheckBound, rightPart, M)

    show_3d_graph(X, Y, Z, center=(x0, y0))


if __name__ == "__main__":
    main()