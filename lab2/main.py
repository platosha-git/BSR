from var2 import difference_model
from data_former import *
from plot_drawer import show_3d_graph
import time

Nx, Ny = 50, 50
Cx, Cy = 46, 20

def main():
	start = time.time()

	X, Y, Z = difference_model(Nx, Ny, Cx, Cy)

	a, b, cc, d = get_data_from_file()
	in_area = form_in_area(a, b, cc, d, Nx, Ny)
	define_in_area(Z, in_area, X[1] - X[0], Y[1] - Y[0])

	end = time.time() - start
	print(end)

	show_3d_graph(X, Y, Z, center=(X[Cx-1], Y[Cy-1], Z[Cx-1][Cy-1]))

if __name__ == "__main__":
	main()