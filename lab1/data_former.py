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