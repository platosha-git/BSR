from plot_drawer import show_3d_graph
from stat_solver import *
from particle import Particle
from random import random

import ctypes
from multiprocessing.sharedctypes import RawArray
import numpy as np
from time import time
from joblib import Parallel, delayed, dump, load

# Make data.

# внешний прямоугольник
a = 60
b = 40

# внутренный прямоугольник
c = 40
d = 20

shiftX = (a - c) / 2 # в центре отверстие
shiftY = (b - d) / 2 

u0 = 300

# f
x0 = a / 2 # shiftX / 2 #
y0 = shiftY / 2 # b / 2 #
alpha = 0.1
k = 0.0134 # 1 #

def f(x, y):
    return 10 * np.exp(-alpha * ((x - x0)**2 + (y-y0)**2))

def rightPart(x, y):
    return f(x, y) / k

X = np.linspace(0, a, 50)
Y = np.linspace(0, b, 50)

hx = X[1] - X[0]
hy = Y[1] - Y[0]
hs = hx * hy

lenX = len(X)
lenY = len(Y)

def f_by_index(i, j):
    return f(X[i], Y[j])

def create_shared_memory_nparray(data):
    d_size = np.dtype(NP_DATA_TYPE).itemsize * np.prod(ARRAY_SHAPE)

    shm = shared_memory.SharedMemory(create=True, size=d_size, name=NP_SHARED_NAME)
    # numpy array on shared memory buffer
    dst = np.ndarray(shape=ARRAY_SHAPE, dtype=NP_DATA_TYPE, buffer=shm.buf)
    dst[:] = data[:]
    print(f'NP SIZE: {(dst.nbytes / 1024) / 1024}')
    return shm


def release_shared(name):
    shm = shared_memory.SharedMemory(name=name)
    shm.close()
    shm.unlink()  # Free and release the shared memory block

raw_z = RawArray(ctypes.c_double, lenX * lenY)
Z = np.frombuffer(raw_z, dtype=np.dtype(raw_z)).reshape(lenX, lenY)
print(raw_z)

Z.fill(-1)
# Z = np.full((len(X), len(Y)), -1, dtype=np.double)

for i in range(len(X)):
    Z[i][0] = u0
    Z[i][-1] = u0

for j in range(len(Y)):
    Z[0][j] = u0
    Z[-1][j] = u0

iShiftX = int(shiftX / hx)
jShiftY = int(shiftY / hy)

indC = int((lenX / a) * c)
indD = int((lenY / b) * d)

print(shiftY, shiftY + d)
print(Y[jShiftY], Y[jShiftY + indD])
print()

print(shiftX, shiftX + c)
print(X[iShiftX], X[iShiftX + indC])
print()

for i in range(iShiftX, iShiftX + indC + 1):
    for j in range(jShiftY, jShiftY + indD):
        Z[i][j] = np.nan

for j in range(jShiftY, jShiftY + indD + 1):
    Z[iShiftX][j] = u0
    Z[iShiftX + indC][j] = u0

for i in range(iShiftX, iShiftX + indC + 1):
    Z[i][jShiftY] = u0
    Z[i][jShiftY + indD] = u0


# краевое 2ого рода
# -k du/dx = F0

# краевое 3его рода
# -k du/dx = a1 (u - u0)

F0 = 30
a0 = 0.1
def CheckBound(p: Particle):
    # левая граница внешнего прямоугольника
    if p.i <= 0:
        # return p.partialReflection(u0, f_by_index, a0, k, hx, di = +1)
        return p.reflection(F0, k, f_by_index, hx, di = +1) # краевое 2ого рода
        # return p.absorption(u0) # краевое 1ого рода

    # правая граница внешнего прямоугольника
    if p.i >= lenX - 1:
        # return p.partialReflection(u0, f_by_index, a0, k, hx, di = -1)
        return p.reflection(F0, k, f_by_index, hx, di = -1)
        # return p.absorption(u0)

    # нижняя граница внешнего прямоугольника
    if p.j <= 0:
        # return p.partialReflection(u0, f_by_index, a0, k, hy, dj = +1)
        return p.reflection(-F0, k, f_by_index, hy, dj = +1)
        # return p.absorption(u0)

    # верхняя граница внешнего прямоугольника
    if p.j >= lenY - 1:
        # return p.partialReflection(u0, f_by_index, a0, k, hy, dj = -1)
        return p.reflection(-F0, k, f_by_index, hy, dj = -1)
        # return p.absorption(u0)

    # во внутреннем прямоугольнике
    if (jShiftY <= p.j and p.j <= jShiftY + indD) and (iShiftX <= p.i and p.i <= iShiftX + indC):
        # левая граница внутреннего прямоугольника
        if (iShiftX == p.i):
            return p.partialReflection(u0, f_by_index, a0, k, hx, di = -1)
            # return p.absorption(u0)

        # правая граница внутреннего прямоугольника
        if (p.i == iShiftX + indC):
            return p.partialReflection(u0, f_by_index, a0, k, hx, di = +1)
            # return p.absorption(u0)

        # нижняя граница внутреннего прямоугольника
        if (jShiftY == p.j):
            return p.partialReflection(u0, f_by_index, a0, k, hy, dj = -1)
            # return p.absorption(u0 - 1)

        # верхняя граница внутреннего прямоугольника
        if (p.j == jShiftY + indD):
            return p.partialReflection(u0, f_by_index, a0, k, hy, dj = +1)
            # return p.reflection(F0, k, f_by_index, hy, dj = +1)
            # return p.absorption(u0)

        # отверстие внутреннего прямоугольника
        p.absorption(np.nan)

# run stochastic method
M = 100

t = time()
# stochastic_solve_point(1, 1, X, Y, Z, CheckBound, rightPart, M)
for i in range(lenX):
    for j in range(lenY):
        stochastic_solve_point(i, j, X, Y, Z, CheckBound, rightPart, M)
print(f"seq stochastic method: {time() - t} sec")

# t = time()

# res = Parallel(n_jobs=4)(delayed(stochastic_solve_point)(i, j, X, Y, Z, IsBound, rightPart, M) for j in range(lenY) for i in range(lenX))
# Z = np.asarray(res).reshape(lenX, lenY)

# print(f"par stochastic method: {time() - t} sec")

zmax = Z[0][0]
imax = 0
jmax = 0
for i in range(len(X)):
    for j in range(len(Y)):
        if zmax < Z[i][j]:
            zmax = Z[i][j]
            imax = i
            jmax = j

print(x0, y0)
print(X[imax], Y[jmax])
print("max Z = ", zmax)
test = (X[imax], Y[jmax], zmax)
np.set_printoptions(precision=0)
# print(Z)

X, Y = np.meshgrid(X, Y, indexing='ij')
show_3d_graph(X, Y, Z, center=(x0, y0)) # test=test