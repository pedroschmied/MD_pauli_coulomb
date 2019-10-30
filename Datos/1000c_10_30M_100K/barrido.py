import numpy as np
import matplotlib.pyplot as plt

M = np.loadtxt('barrido_rho_0.15.txt', delimiter = '\t')
T = M[:, 0]
E = M[:, 1]
aceptacion_x = M[:, 2]
aceptacion_p = M[:, 3]

M2 = np.loadtxt('barrido_rho_0.14.txt', delimiter = '\t')
T2 = M2[:, 0]
E2 = M2[:, 1]
aceptacion_x2 = M2[:, 2]
aceptacion_p2 = M2[:, 3]
'''
M3 = np.loadtxt('corrida_T_4.00_rho_0.14.txt', delimiter = '\t')
T3 = M3[:, 0]
E3 = M3[:, 1]
#aceptacion_x3 = M3[:, 2]
#aceptacion_p3 = M3[:, 3]
'''
plt.figure()

plt.grid()
plt.ylabel('E (MeV)', fontsize = 13)
plt.xlabel('T (MeV)', fontsize = 13)
plt.title('Curva calórica', fontsize = 15)

plt.plot(T, E, '.', label = 'rho = 0.15')
plt.plot(T2, E2, '.', label = 'rho = 0.14')
#plt.plot(T3, E3, '.', label = 'rho = 0.14')
#plt.plot(x, np.sqrt(x) * ap + bp, label = 'ajuste dp')
#plt.plot(T, dx, '.', label = 'rho = 0.10')
plt.legend()
plt.show()

