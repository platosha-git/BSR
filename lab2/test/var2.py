from math import exp
from data_former import *
from plot_drawer import show_3d_graph

sigma = 5.669e-8;
eps = 1e-5;

Nx = 50
Ny = 50
t_end = 36000.0
L = 0.6
H = 0.4
lamda = 0.16
ro = 1190.0
c = 1900.0
kapa1 = 50.0
kapa2 = 35.0
Te1 = 200.0
Te2 = 200.0
eps1 = 0.8
T0 = 300.0
q = 10000.0
alpha = 0.1

def sqr(x):
	return x**2

def main():
	#определяем расчетные шаги сетки по пространственным координатам
	hx = L / (Nx - 1);
	hy = H / (Ny - 1);

	N1x = 46
	N1y = 20

	#N1x = 4
	#N1y = 1
	
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
			alfa[1] = 2.0 * tau * lamda / (2.0 * tau * (lamda + kapa1 * hx) + ro * c * sqr(hx));

			#вычисляем поле температуры, вследствие наличия нелинейности в левом граничном условии
			while True:
				#определяем beta начальный прогоночный коэффициент на основе левого граничного условия 
				#начинаем итерационный цикл по левому граничному условию
		  
				d = T[1][j];
				beta[1] = (ro * c * sqr(hx) * Tn[1][j] + \
								2.0 * tau * kapa1 * hx * Te1 + \
								2.0 * tau * eps1 * sigma * hx * (sqr(sqr(Te1))-sqr(sqr(d)))) / \
							(2.0 * tau * (lamda + kapa1 * hx) + ro * c * sqr(hx));	

				#цикл для определения прогоночных коэффициентов
				for i in range(2, Nx):
					
					#ai, bi, ci, fi – коэффициенты канонического представления СЛАУ с трехдиагональной матрицей
					ai = lamda / sqr(hx);
					bi = 2.0 * lamda / sqr(hx) + ro * c / tau;
					ci = lamda / sqr(hx);
					fi = -ro * c * Tn[i][j] / tau;

					if i == N1x:
						fi = fi - q * exp(-alpha * (sqr(i - N1x)));

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

				if j == N1y:
					fi = fi - q * exp(-alpha * (sqr(j - N1y)));

				#alfa[j], beta[j] – прогоночные коэффициенты
				alfa[j] = ai / (bi - ci * alfa[j-1]);
				beta[j] = (ci * beta[j-1] - fi) / (bi - ci * alfa[j-1]);			


			#запоминаем значение температуры на правой границе с промежуточного (n+1/2) временного слоя
			d = T[i][Ny];

			while True:
				d1 = T[i][Ny];

				#определяем значение температуры на правой границе на основе правого граничного условия
				T[i][Ny] = (ro * c * sqr(hy) * d + 2.0 * tau * \
								(lamda * beta[Ny-1] + kapa2 * hy * Te2 + eps1 * sigma * hy * (sqr(sqr(Te2))-sqr(sqr(d1))))) / \
							(ro * c * sqr(hy) + 2.0 * tau * (lamda * (1 - alfa[Ny-1]) + kapa2 * hy));
				
				if abs((d1 - T[i][Ny]) / T[i][Ny]) <= eps:
					break

			for j in range(Ny - 1, -1, -1):	
				T[i][j] = alfa[j] * T[i][j+1] + beta[j];


	print('Длина пластины L = ', L);
	print('Толщина пластины H = ', H);
	print('Число узлов по пространственной координате x в пластине Nx = ',Nx);
	print('Число узлов по пространственной координате y в пластине Ny = ',Ny);
	print('\n')
	print('Начальная температура T0 = ', T0);
	print('Температура внешней среды Te1 = ', Te1);
	print('\n');
	print('Результат получен с шагом по координате x hx = ', hx);
	print('Результат получен с шагом по координате y hy = ', hy);
	print('\n');
	print('Результат получен с шагом по времени tau = ', tau);
	print('Температурное поле в момент времени t = ', t_end);

	
	X = [0.0 for i in range(Nx)]
	for i in range(1, Nx+1):
		X[i-1] = hx * (i-1)

	Y = [0.0 for j in range(Ny)]
	for j in range(1, Ny+1):
		Y[j-1] = hy * (j-1)

	Z = [[0.0] * Ny for i in range(Nx)]
	for i in range(1, Nx+1):
		for j in range(1, Ny+1):
			Z[i-1][j-1] = T[i][j]

	a, b, cc, d = get_data_from_file()
	in_area = form_in_area(a, b, cc, d, Nx, Ny)
	define_in_area(Z, in_area, X[1] - X[0], Y[1] - Y[0])

	show_3d_graph(X, Y, Z, center=(X[46], Y[20], Z[46][20]))

if __name__ == "__main__":
	main()