from plot_drawer import show_3d_graph
from stat_solver import *
from data_former import *
from particle import Particle
from random import random
import numpy as np
from time import time
from joblib import Parallel, delayed, dump, load


k = 0.0134 # 0.001 # 1 #

def f(x, y, x0, y0):
    alpha = 0.1
    return 10 * np.exp(-alpha * ((x - x0)**2 + (y-y0)**2))

def f_by_index(i, j, X, Y, x0, y0):
    return f(X[i], Y[j], x0, y0)

def rightPart(x, y, x0, y0):
    return f(x, y, x0, y0) / k


F0 = 30
a0 = 0.1
def CheckBound(p: Particle, X, Y, x0, y0, cX, cY, indC, indD, u0):
    lenX = len(X)
    lenY = len(Y)

    hX = X[1] - X[0]
    hY = Y[1] - Y[0]

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
    
    cX = (a - c) / 2
    cY = (b - d) / 2 

    x0 = a / 2 
    # y0 = b / 2
    
    # x0 = cX / 2
    y0 = cY / 2

    indC = int((lenX / a) * c)
    indD = int((lenY / b) * d)

    #print("x=", X)
    #print("y=", Y)
    #print("z=", Z)

    #print("cX = ", cX)
    #print("cY = ", cY)

    #hS = hX * hY

    #print("hX = ", hX)
    #print("hY = ", hY)

    #print("icX = ", icX)
    #print("jcY = ", jcY)

    #print(cY, cY + d)
    #print(Y[jShiftY], Y[jShiftY + indD])
    #print()

    #print(shiftX, shiftX + c)
    #print(X[iShiftX], X[iShiftX + indC])
    #print()

    # run stochastic method
    M = 10
    t = time()
    for i in range(lenX):
        for j in range(lenY):
            stochastic_solve_point(i, j, X, Y, Z, x0, y0, cX, cY, indC, indD, u0, CheckBound, rightPart, M)
    print(f"seq stochastic method: {time() - t} sec")




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



if __name__ == "__main__":
    main()

'''
# f

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

# Z = np.full((len(X), len(Y)), -1, dtype=np.double)



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
'''