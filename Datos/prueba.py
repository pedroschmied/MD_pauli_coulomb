import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

M = np.loadtxt('barrido_rho_0.16.txt', delimiter = '\t')
pasos = M[:, 0]
E = M[:, 1]

#M1 = np.loadtxt('test_tablas_V_LJ.txt', delimiter = '\t')
#pasos1 = M1[:, 0]
#E1 = M1[:, 1]

plt.figure()
plt.grid()
#plt.ylim((-5.4,2.0))
plt.plot(pasos, E, '.', label = 'hola')
#plt.plot(pasos1, E1, '.', label = 'Pedro')

plt.show()
plt.legend()
