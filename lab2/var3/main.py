from plot_drawer import show_3d_graph
from data_former import *

X_filename = './X.txt'
Y_filename = './Y.txt'
Z_filename = './Z.txt'

X_size = 50
Y_size = 50

def main():
	X = read_vector(X_filename)
	Y = read_vector(Y_filename)
	Z = read_matrix(Z_filename, X_size, Y_size)

	a, b, c, d = get_data_from_file()
	in_area = form_in_area(a, b, c, d, X_size, Y_size)
	define_in_area(Z, in_area, X[1] - X[0], Y[1] - Y[0])

	show_3d_graph(X, Y, Z, center=(X[46], Y[20], Z[46][20]))


if __name__ == "__main__":
	main()
	