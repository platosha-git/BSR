from multiprocessing.sharedctypes import RawArray
import numpy as np
import ctypes
import json

filename = 'data.json'

def get_data_from_file():
	with open(filename, 'r') as f:
		text = json.load(f)
		
		outer = text["Outer area"]
		a = outer["a"]
		b = outer["b"]

		inner = text["Inner area"]
		c = inner["c"]
		d = inner["d"]

	return a, b, c, d


def form_plane(a, b):
	X = np.linspace(0, a)
	Y = np.linspace(0, b)
	
	lenX = len(X)
	lenY = len(Y)

	buff = RawArray(ctypes.c_double, lenX * lenY)
	Z = np.frombuffer(buff, dtype=np.dtype(buff)).reshape(lenX, lenY)
	Z.fill(-1)

	return X, Y, Z


def form_circuit(Z, u0):
	Z[0].fill(u0)
	Z[-1].fill(u0)

	Z[:, 0].fill(u0)
	Z[:, -1].fill(u0)



inner_area = {
	'border': [0, 0],
	'size': [0, 0]
}

# Отверстие расположено посередине
def form_in_area(a, b, c, d, lenX, lenY):
	inner_area['border'][0] = (a - c) / 2
	inner_area['border'][1] = (b - d) / 2

	inner_area['size'][0] = int((lenX / a) * c)
	inner_area['size'][1] = int((lenY / b) * d)

	return inner_area


def define_in_area(Z, in_area, hX, hY):
	lX = int(in_area['border'][0] / hX)
	lY = int(in_area['border'][1] / hY)

	rX = lX + inner_area['size'][0] + 1
	rY = lY + inner_area['size'][1]

	for i in range(lX, rX):
		for j in range(lY, rY):
			Z[i][j] = np.nan


def define_output_area(Z, in_area, hX, hY, u0):
	lX = int(in_area['border'][0] / hX)
	lY = int(in_area['border'][1] / hY)

	rX = lX + inner_area['size'][0]
	rY = lY + inner_area['size'][1]
	
	for j in range(lY, rY + 1):
		Z[lX][j] = u0
		Z[rX][j] = u0

	for i in range(lX, rX + 1):
		Z[i][lY] = u0
		Z[i][rY] = u0