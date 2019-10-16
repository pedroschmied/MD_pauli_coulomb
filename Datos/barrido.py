import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

M = np.loadtxt('barrido_rho_0.16.txt', delimiter = '\t')
T = M[:, 0]
E = M[:, 1]
Cv = M[:, 2]
aceptacion_x = M[:, 3]
aceptacion_p = M[:, 4]
aceptacion = M[:, 5]

M2 = np.loadtxt('barrido_rho_0.10.txt', delimiter = '\t')
T2 = M2[:, 0]
E2 = M2[:, 1]
Cv2 = M2[:, 2]
aceptacion_x2 = M2[:, 3]
aceptacion_p2 = M2[:, 4]
aceptacion2 = M2[:, 5]

#M1 = np.loadtxt('test_tablas_V_LJ.txt', delimiter = '\t')
#pasos1 = M1[:, 0]
#E1 = M1[:, 1]

plt.figure()
plt.grid()
plt.ylabel('E (MeV)', fontsize = 13)
plt.xlabel('T (MeV)', fontsize = 13)
plt.title('Curva calórica', fontsize = 15)

plt.plot(T, aceptacion, '.', label = 'rho = 0.16')
plt.plot(T2, aceptacion2, '.', label = 'rho = 0.10')
plt.legend()
plt.show()

