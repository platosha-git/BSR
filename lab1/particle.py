from random import random
import numpy as np
from typing import Callable

class Particle:
    def __init__(self, i, j):
        self.i = i
        self.j = j
        self.isAlive = True
        self.u = np.nan

        self.A = 0
        self.n = 0

    def addA(self, a):
        self.A += a
        self.n += 1

    def getA(self):
        if self.n > 0:
            return self.A / self.n
        return self.A

    def next(self):
        r = random()
        if r < 0.25:
            self.i += 1
        elif r < 0.5:
            self.j += 1
        elif r < 0.75:
            self.i -= 1
        else:
            self.j -= 1
    
    # поглощение на границе со значением u0 (краевая задача 1ого рода)
    def absorption(self, u0):
        self.isAlive = False
        self.u = u0

    # отражение на границе (краевая задача 2ого рода -- задан поток тепла F0)
    def reflection(self, F0, k, Q: Callable[[int, int], float], h, di=0, dj=0):
        self.i += di
        self.j += dj
        self.addA((Q(self.i, self.j) * h / 2 + F0) * h / k)

    # вероятностное отражение на границе (краевая задача 3ого рода)
    # u0 - начальная температура
    # a0 - коэф теплоотдачи
    def partialReflection(self, u0, Q: Callable[[int, int], float], a0, k, h, di=0, dj=0):
        r = random()
        
        # вероятность отражения
        p0 = 1 / (1 + a0 / k * h)
        if r < p0:
            self.i += di
            self.j += dj
            self.addA(p0 * Q(self.i, self.j) * h * h / 2 / k)
        else:
            # иначе поглощение
            self.absorption(u0)