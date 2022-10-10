import numpy as np
import numpy.typing as npt
from typing import Callable
from particle import Particle
from main import f, f_by_index

def stochastic_solve_point(i, j, X, Y, Z, heat_point, in_area, u0, checkBound, f, M = 100):
    u = 0
    for _ in range(M):
        p = Particle(i, j)
        checkBound(p, X, Y, heat_point, in_area, u0)
        
        while p.isAlive:
            p.addA(f(X[p.i], Y[p.j], heat_point[0], heat_point[1]))
            p.next()
            
            checkBound(p, X, Y, heat_point, in_area, u0)

        ak = p.getA()
        u += ak + p.u

    u = u / M
    Z[i][j] = u
    
    return u


def rightPart(x, y, x0, y0):
    return f(x, y, x0, y0) / k

# краевое 2ого рода
# -k du/dx = F0

# краевое 3его рода
# -k du/dx = a1 (u - u0)

k = 0.01
F0 = 30
a0 = 0.1

def CheckBound(p: Particle, X, Y, heat_point, in_area, u0):
    lenX = len(X)
    lenY = len(Y)

    hX = X[1] - X[0]
    hY = Y[1] - Y[0]

    lX = int(in_area['border'][0] / hX)
    lY = int(in_area['border'][1] / hY)

    rX = lX + in_area['size'][0]
    rY = lY + in_area['size'][1]


    # Outer area: left, right, down, up
    if p.i <= 0:
        return p.reflection(-F0, k, f_by_index, X, Y, heat_point, hX, di = +1)

    if p.i >= lenX - 1:
        return p.reflection(F0, k, f_by_index, X, Y, heat_point, hX, di = -1)

    if p.j <= 0:
        return p.reflection(F0, k, f_by_index, X, Y, heat_point, hY, dj = +1)

    if p.j >= lenY - 1:
        return p.reflection(F0, k, f_by_index, X, Y, heat_point, hY, dj = -1)


    # Inner area: left, right, down, up
    if (lY <= p.j and p.j <= rY) and (lX <= p.i and p.i <= rX):
        if (lX == p.i):
            return p.partialReflection(u0, f_by_index, X, Y, heat_point, a0, k, hX, di = -1)

        if (p.i == rX):
            return p.partialReflection(u0, f_by_index, X, Y, heat_point, a0, k, hX, di = +1)

        if (lY == p.j):
            return p.partialReflection(u0, f_by_index, X, Y, heat_point, a0, k, hY, dj = -1)

        if (p.j == rY):
            return p.partialReflection(u0, f_by_index, X, Y, heat_point, a0, k, hY, dj = +1)

        p.absorption(np.nan)
