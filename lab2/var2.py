from math import exp

sigma = 5.669e-8;
eps = 1e-5;

t_end = 36000.0
L = 0.6
H = 0.4
lamda = 0.16
ro = 1190.0
c = 1900.0
kapa = 50.0
Te = 200.0
eps1 = 0.8
T0 = 300.0

def output_result(Nx, Ny, hx, hy, tau):
	print('Длина пластины L = ', L);
	print('Толщина пластины H = ', H);
	print('Число узлов по пространственной координате x в пластине Nx = ', Nx);
	print('Число узлов по пространственной координате y в пластине Ny = ', Ny);
	print('\n')
	print('Начальная температура T0 = ', T0);
	print('Температура внешней среды Te = ', Te);
	print('\n');
	print('Результат получен с шагом по координате x hx = ', hx);
	print('Результат получен с шагом по координате y hy = ', hy);
	print('\n');
	print('Результат получен с шагом по времени tau = ', tau);
	print('Температурное поле в момент времени t = ', t_end);


def form_result(hx, hy, Nx, Ny, T):
	X = [0.0 for i in range(Nx)]
	for i in range(Nx):
		X[i] = hx * i

	Y = [0.0 for j in range(Ny)]
	for j in range(Ny):
		Y[j] = hy * j

	Z = [[0.0] * Ny for i in range(Nx)]
	for i in range(Nx):
		for j in range(Ny):
			Z[i][j] = T[i+1][j+1]

	return X, Y, Z


def f(x, x0, y, y0):
	A = 100000.0
	alpha = 2
	return A * exp(-alpha * (sqr(x - x0) + sqr(y - y0)));


def sqr(x):
	return x**2


def difference_model(Nx, Ny, Cx, Cy):
	#определяем расчетные шаги сетки по пространственным координатам
	hx = L / (Nx - 1);
	hy = H / (Ny - 1);

	N1x = Cx
	N1y = Cy
	
	#определяем расчетный шаг сетки по времени
	tau = t_end / 1000.0;

	#определяем коэффициент температуропроводности
	a = lamda/(ro*c)
	
	T = [[T0] * (Ny+1) for i in range(Nx+1)]
	Tn = [[0.0] * (Ny+1) for i in range(Nx+1)]

	alfa = [0.0 for i in range(Nx+1)]
	beta = [0.0 for i in range(Ny+1)]

	#проводим интегрирование нестационарного уравнения теплопроводности
	time = 0.0;
	while time < t_end:
		time = time + tau;

		#запоминаем поле температуры на n-ом временном слое
		for i in range(1, Nx+1):
			for j in range(1, Ny+1):
				Tn[i][j] = T[i][j]

		#решаем СЛАУ в направлении оси Ох для определения поля температуры на промежуточном (n+1/2) временном слое
		for j in range(1, Ny+1):
			alfa[1] = 2.0 * tau * lamda / (2.0 * tau * (lamda + kapa * hx) + ro * c * sqr(hx));

			#вычисляем поле температуры, вследствие наличия нелинейности в левом граничном условии
			while True:
				#определяем beta начальный прогоночный коэффициент на основе левого граничного условия 
				#начинаем итерационный цикл по левому граничному условию
		  
				d = T[1][j];
				beta[1] = (ro * c * sqr(hx) * Tn[1][j] + \
								2.0 * tau * kapa * hx * Te + \
								2.0 * tau * eps1 * sigma * hx * (sqr(sqr(Te))-sqr(sqr(d)))) / \
							(2.0 * tau * (lamda + kapa * hx) + ro * c * sqr(hx));	

				#цикл для определения прогоночных коэффициентов
				for i in range(2, Nx):
					
					#ai, bi, ci, fi – коэффициенты канонического представления СЛАУ с трехдиагональной матрицей
					ai = lamda / sqr(hx);
					bi = 2.0 * lamda / sqr(hx) + ro * c / tau;
					ci = lamda / sqr(hx);
					fi = -ro * c * Tn[i][j] / tau;

					if i == N1x and j == N1y:
						fi = fi - f(i, N1x, j, N1y)

					#alfa[i], beta[i] – прогоночные коэффициенты
					alfa[i] = ai / (bi - ci * alfa[i - 1]);
					beta[i] = (ci * beta[i-1] - fi) / (bi - ci * alfa[i-1]);
				
				#определяем значение температуры на правой границе
				T[Nx][j] = (ro * c * sqr(hx) * Tn[Nx][j] + 2.0 * tau * lamda * beta[Nx - 1]) / \
							(ro * c * sqr(hx) + 2.0 * tau * lamda * (1 - alfa[Nx - 1]));
				
				#определяем неизвестное поле температуры на промежуточном (n+1/2) временном слое
				for i in range(Nx - 1, 0, -1):
					T[i][j] = alfa[i] * T[i+1][j] + beta[i];
			
				if abs((d - T[1][j]) / T[1][j]) <= eps:
					break
		
		
		#поле температуры на промежуточном (n+1/2) временном слое определили
		#решаем СЛАУ в направлении оси Оу для определения поля температуры на целом (n+1) временном слое
		for i in range(1, Nx+1):

			#определяем начальные прогоночные коэффициенты на основе нижнего граничного условия
			alfa[1] = 2.0 * tau * lamda / (2.0 * tau * lamda + ro * c * sqr(hy));
			beta[1] = ro * c * sqr(hy) * T[i][1] / (2.0 * tau * lamda + ro * c * sqr(hy));

			#цикл для определения прогоночных коэффициентов
			for j in range(2, Ny):

				#ai, bi, ci, fi – коэффициенты канонического представления СЛАУ с трехдиагональной матрицей
				ai = lamda / sqr(hy);
				bi = 2.0 * lamda / sqr(hy) + ro * c / tau;
				ci = lamda / sqr(hy);
				fi = -ro * c * T[i][j] / tau;

				if i == N1x and j == N1y:
					fi = fi - f(i, N1x, j, N1y)

				#alfa[j], beta[j] – прогоночные коэффициенты
				alfa[j] = ai / (bi - ci * alfa[j-1]);
				beta[j] = (ci * beta[j-1] - fi) / (bi - ci * alfa[j-1]);			


			#запоминаем значение температуры на правой границе с промежуточного (n+1/2) временного слоя
			d = T[i][Ny];

			while True:
				d1 = T[i][Ny];

				#определяем значение температуры на правой границе на основе правого граничного условия
				T[i][Ny] = (ro * c * sqr(hy) * d + 2.0 * tau * \
								(lamda * beta[Ny-1] + kapa * hy * Te + eps1 * sigma * hy * (sqr(sqr(Te))-sqr(sqr(d1))))) / \
							(ro * c * sqr(hy) + 2.0 * tau * (lamda * (1 - alfa[Ny-1]) + kapa * hy));
				
				if abs((d1 - T[i][Ny]) / T[i][Ny]) <= eps:
					break

			for j in range(Ny - 1, -1, -1):	
				T[i][j] = alfa[j] * T[i][j+1] + beta[j];


	output_result(Nx, Ny, hx, hy, tau)
	X, Y, Z = form_result(hx, hy, Nx, Ny, T)

	return X, Y, Z