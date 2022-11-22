from plot_drawer import show_3d_graph

sigma = 5.669e-8;
eps = 1e-5;

Nx = 5
Ny = 5
t_end = 36000
L = 0.6
H = 0.4
lamda = 0.16
ro = 1190.0
c = 1900.0
kapa1 = 50.0
kapa2 = 35.0
Te1 = 20.0
Te2 = 35.0
eps1 = 0.8
T0 = 30.0

def f(x, y, x0, y0):
    alpha = 2
    res = 10 * np.exp(-alpha * ((x - x0)**2 + (y-y0)**2))
    return res

def sqr(x):
	return x**2

def main():
	hx = L / (Nx - 1);
	hy = H / (Ny - 1);
	
	tau = t_end / 1000.0;

	T = [[T0] * (Ny+1) for i in range(Nx+1)]
	Tn = [[0] * (Ny+1) for i in range(Nx+1)]

	alfa = [0.0 for i in range(Nx+1)]
	beta = [0.0 for i in range(Ny+1)]

	time = 0;
	while time < tau:
		time = time + tau;

		for i in range(1, Nx+1):
			for j in range(1, Ny+1):
				Tn[i][j] = T[i][j]

		for j in range(1, Ny+1):
			alfa[1] = 2.0 * tau * lamda / (2.0 * tau * (lamda + kapa1 * hx) + ro * c * pow(hx, 2));

			while True:
				d = T[1][j];
				beta[1] = (ro * c * sqr(hx) * Tn[1][j] + 2.0 * tau * kapa1 * hx * Te1 + \
							2.0 * tau * eps1 * sigma * hx * (sqr(sqr(Te1)) - sqr(sqr(d)))) /\
							(2.0 * tau * (lamda + kapa1 * hx) + ro * c * sqr(hx));	

				for i in range(2, Nx):
					ai = lamda / pow(hx, 2);
					bi = 2.0 * lamda / pow(hx, 2) + ro * c / tau;
					ci = lamda / pow(hx, 2);
					fi = -ro * c * Tn[i][j] / tau;

					if (j == N1y):
						fi = -ro * c * T[i] / tau - f(i, j, x0, y0);

					alfa[i] = ai / (bi - ci * alfa[i - 1]);
					beta[i] = (ci*beta[i-1]-fi)/(bi-ci*alfa[i-1]);
				
				T[Nx][j] = (ro * c * pow(hx, 2) * Tn[Nx][j] + 2.0 * tau * lamda * beta[Nx - 1]) / \
							(ro * c * pow(hx, 2) + 2.0 * tau * lamda * (1 - alfa[Nx - 1]));
				
				for i in range(Nx - 1, 0, -1):
					T[i][j] = alfa[i] * T[i+1][j] + beta[i];
			
				if abs((d - T[1][j]) / T[1][j]) <= eps:
					break
		
		for i in range(Nx+1):
			alfa[1] = 2.0 * tau * lamda / (2.0 * tau * lamda + ro * c * sqr(hy));
			beta[1] = ro*c*sqr(hy)*T[i][1]/(2.0*tau*lamda+ro*c*sqr(hy));

			for j in range(2, Ny):
				ai = lamda / sqr(hy);
				bi = 2.0 * lamda / sqr(hy) + ro * c / tau;
				ci = lamda / sqr(hy);
				fi = -ro * c * T[i][j] / tau;

				if (j == N1y):
					fi = -ro * c * T[i] / tau - f(i, j, x0, y0);

				alfa[j] = ai / (bi - ci * alfa[j-1]);
				beta[j] = (ci * beta[j-1] - fi) / (bi - ci * alfa[j-1]);			
			
			d = T[i][Ny];

			while True:
				d1 = T[i][Ny];
				T[i][Ny] = (ro * c * sqr(hy) * d + 2.0 * tau * (lamda * beta[Ny-1] + kapa2 * hy * Te2 + \
							eps1 * sigma * hy * (sqr(sqr(Te2)) - sqr(sqr(d1))))) / \
							(ro * c * sqr(hy) + 2.0 * tau * (lamda * (1 - alfa[Ny-1])+kapa2*hy));
				
				if abs((d1 - T[i][Ny]) / T[i][Ny]) <= eps:
					break

			for j in range(Ny - 1, 0, -1):	
				T[i][j] = alfa[j] * T[i][j+1] + beta[j];

	print('Длина пластины L = ',L);
	print('Толщина пластины H = ',H);
	print('Число узлов по пространственной координате x в пластине Nx = ',Nx);
	print('Число узлов по пространственной координате y в пластине Ny = ',Ny);
	print('Начальная температура T0 = ',T0);
	print('Температура внешней среды Te1 = ',Te1);
	print('Приведенная степень черноты eps1 = ',eps1);
	print('Результат получен с шагом по координате x hx = ',hx);
	print('Результат получен с шагом по координате y hy = ',hy);
	print('Результат получен с шагом по времени tau = ',tau);
	print('Температурное поле в момент времени t = ',t_end);

	X = [0 for i in range(Nx)]
	Y = [0 for j in range(Ny)]
	Z = [[0] * Ny for i in range(Nx)]

	for i in range(Nx):
		X[i] = hx * i
	
	for j in range(Ny):
		Y[j] = hy * j


	for i in range(Nx):
		for j in range(Ny):
			print(' ', hx*(i),' ', hy*(j),' ',T[i][j]);
	# 		Z[i][j] = T[i][j] * 100

	#show_3d_graph(X, Y, T, center=(heat_point[0], heat_point[1]))
	#show_3d_graph(X, Y, Z)

if __name__ == "__main__":
	main()