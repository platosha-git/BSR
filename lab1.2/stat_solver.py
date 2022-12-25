from random import random
import numpy as np
import numpy.typing as npt
from typing import Callable
from particle import Particle

def stochastic_solve_point(i: int, j: int, \
        X: npt.ArrayLike, Y: npt.ArrayLike, Z: npt.NDArray[float], \
        checkBound: Callable[[Particle], bool], f: Callable[[float, float], float], M = 100):
    u = 0
    for _ in range(M):
        p = Particle(i, j)
        checkBound(p)
        while p.isAlive:
            p.addA(f(X[p.i], Y[p.j]))

            p.next()
            checkBound(p)

        ak = p.getA()
        u += ak + p.u

    u = u / M
    Z[i][j] = u
    return u

def stochastic_solve(X: npt.ArrayLike, Y: npt.ArrayLike, Z: npt.NDArray[float], \
        isBound: Callable[[int, int], bool], f: Callable[[float, float], float], M = 100):
    lenX = len(X)
    lenY = len(Y)

    for i in range(lenX):
        for j in range(lenY):
            stochastic_solve_point(i, j, X, Y, Z, isBound, f, M)
    return Z