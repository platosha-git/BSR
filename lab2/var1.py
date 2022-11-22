import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

N = 100
t_end = 10
L = 0.3
lamda = 419.0
ro = 10500.0
c = 200.0
Te = 60.0
kapa = 50.0
T0 = 300
q = 100000.0

def sqr(x):
	return x**2

def main():
	#определяем узел, в котором расположен источник
	N1 = N / 4

	#определяем расчетный шаг сетки по пространственной координате
	h = L / (N - 1)

	#определяем коэффициент температуропроводности
	a = lamda / (ro * c)

	#определяем расчетный шаг сетки по времени
	tau = t_end / 100.0

	T = [T0 for i in range(N+1)]
	alfa = [0.0 for i in range(N+1)]
	beta = [0.0 for i in range(N+1)]

	time = 0.0
	while time < t_end:
		time = time + tau

		#определяем начальные прогоночные коэффициенты на основе левого граничного условия
		alfa[1] = 2.0 * a * tau * lamda / (lamda * sqr(h) +\
					2.0 * a * tau * (lamda + kapa *h));
		beta[1] = (lamda * sqr(h) * T[1] + 2.0 * a * tau * kapa * h * Te) / \
					(lamda * sqr(h) + 2.0 * a * tau * (lamda + kapa * h));

		for i in range(2, N):
			#ai, bi, ci, fi – коэффициенты канонического представления системы уравнений с трехдиагональной матрицей
			ai = lamda / sqr(h);
			bi = 2.0 * lamda / sqr(h) + ro * c / tau;
			ci = lamda / sqr(h);
			fi = -ro * c * T[i] / tau;

			#определяем fi в зависимости от рассматрвиаемой точки пространства
			if i == N1:
				fi = -ro * c * T[i] / tau - h * (i - 1) * q;

			#alfa[i], beta[i] – прогоночные коэффициенты
			alfa[i] = ai / (bi - ci * alfa[i-1]);
			beta[i] = (ci * beta[i-1] - fi) / (bi - ci * alfa[i-1]);

		#определяем значение температуры на правой границе
		T[N] = (lamda * sqr(h) * T[N] + 2.0 * a * tau * (lamda * beta[N-1] + kapa * h * Te)) /\
					(lamda * sqr(h) + 2.0 * a * tau * (lamda * (1 - alfa[N-1]) + kapa * h));

		for i in range(N - 1, -1, -1):
			T[i] = alfa[i] * T[i + 1] + beta[i]


	print('Толщина пластины L = ', L)
	print('Число узлов по пространственной координате в пластине N = ', N)
	print('Внутренний источник находятся в точке = ', N1)
	print('Внутренний источник находятся в точке = ', (N1-1)*h)
	print('\nНачальная температура T0 = ', T0)
	print('Температура окружающей среды Te = ', Te)
	print('Составляющая мощности внутренних источников тепла q = ', q)
	print('\nРезультат получен с шагом по координате h = ', h)
	print('Результат получен с шагом по времени tau = ', tau)
	print('Температурное поле в момент времени t = ', t_end)

	x = [0 for i in range(N+1)]
	for i in range(1, N+1):
		x[i] = h * (i-1)
		#print(' ',x[i],' ', T[i]);


	fig, ax = plt.subplots()
	ax.plot(x, T)
	
	plt.show()


if __name__ == "__main__":
	main()
