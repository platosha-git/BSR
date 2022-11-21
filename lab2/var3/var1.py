from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

mf = 500;
vector = [i + 1 for i in range(mf)]


def main():
	N = 100
	t_end = 1
	L = 0.3
	lamda = 419.0
	ro = 10500.0
	c = 200.0
	Te = 60.0
	kapa = 50.0
	T0 = 300
	q = 100000.0

	#{определяем номера узлов, в которых расположен источник}
	N1 = N / 4

	#{определяем расчетный шаг сетки по пространственной координате}
	h = L / (N - 1)
	a = lamda / (ro * c)

	tau = t_end / 100.0

	T = [T0 for i in range(N)]
	alfa = [0.0 for i in range(N)]
	beta = [0.0 for i in range(N)]

	time = 0.0
	while time < t_end:
		time = time + tau

		#{определяем начальные прогоночные коэффициенты на основе левого граничного условия, используя соотношения (24)}
		alfa[0] = 2.0 * a * tau * lamda / (lamda * sqrt(h) +\
					2.0 * a * tau * (lamda + kapa *h));
		beta[0] = (lamda * sqrt(h) * T[0] + 2.0 * a * tau * kapa * h * Te) / \
					(lamda * sqrt(h) + 2.0 * a * tau * (lamda + kapa * h));

		for i in range(1, N):
			#{ai, bi, ci, fi – коэффициенты канонического представления системы уравнений с трехдиагональной матрицей}
			ai = lamda / sqrt(h);
			bi = 2.0 * lamda / sqrt(h) + ro * c / tau;
			ci = lamda / sqrt(h);
			fi = -ro * c * T[i] / tau;

			#{определяем fi в зависимости от рассматрвиаемой точки пространства}
			if (i == N1):
				fi = -ro * c * T[i] / tau - h * (i-1) * q;

			alfa[i] = ai / (bi - ci * alfa[i-1]);
			beta[i] = (ci * beta[i-1] - fi) / (bi - ci * alfa[i-1]);

		T[N-1] = (lamda * sqrt(h) * T[N-1] + 2.0 * a * tau * (lamda * beta[N-2] + kapa * h * Te)) /\
					(lamda * sqrt(h) + 2.0 * a * tau * (lamda * (1 - alfa[N-2]) + kapa * h));

		for i in range(N - 2, -1, -1):
			T[i] = alfa[i] * T[i + 1] + beta[i]

	'''
	print('Толщина пластины L = ',L)
	print('Число узлов по пространственной координате в пластине N = ',N)
	print('Внутренние источники находятся в точках = ',N1)
	print('Внутренние источники находятся в точках = ',(N1-1)*h)
	print('Коэффициент теплопроводности материала пластины lamda = ',lamda)
	print('Плотность материала пластины ro = ',ro)
	print('Теплоемкость материала пластины с = ',c)
	print('Начальная температура T0 = ',T0)
	print('Коэффициент теплообмена kapa = ',kapa)
	print('Температура окружающей среды Te = ',Te)
	print('Составляющая мощности внутренних источников тепла q = ',q)
	print('Результат получен с шагом по координате h = ',h)
	print('Результат получен с шагом по времени tau = ',tau)
	print('Температурное поле в момент времени t = ',t_end)
	'''

	x = [0 for i in range(N)]
	y = [0 for i in range(N)]
	for i in range(N):
		x[i] = h * i
		y[i] = T[i] * 10000
		print(' ',x[i],' ', T[i]);


	fig, ax = plt.subplots()
	ax.plot(x, y)
	#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
	#  Устанавливаем интервал вспомогательных делений:
	#ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
	
	plt.show()


if __name__ == "__main__":
	main()
