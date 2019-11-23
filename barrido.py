import numpy as np
import matplotlib.pyplot as plt


def barrido(rhof, N, V):
	for i in range (0, N):
		j = rhof - i
		M = np.loadtxt('512_barrido_rho_0.%i.txt' %j, delimiter = '\t')
		T = M[:, 0]
		E = M[:, 1]
		aceptacion_x = M[:, 2]
		aceptacion_p = M[:, 3]
		L = M[:, V]
		plt.plot(T, L, '.', label = 'rho = 0.%i' %j)
	return

plt.figure()

plt.grid()
plt.ylabel('E (MeV)', fontsize = 13)
plt.xlabel('T (MeV)', fontsize = 13)
plt.title('Curva calórica', fontsize = 15)

#barrido(13, 1, 1)
barrido(14, 1, 1)

plt.legend()
plt.show()

