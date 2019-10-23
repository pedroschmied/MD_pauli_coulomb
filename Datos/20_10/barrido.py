import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def ajuste_raiz (logx, logy):

    def ajuste(x,a,b):
        return np.sqrt(x) * a + b

    popt,pcov = curve_fit (ajuste, logx, logy)
    a = popt [0]
    b = popt[1]
#    c = popt[2]
    perr = np.sqrt(np.diag(pcov))

    logymean = np.mean(logy)

    SSTOT=np.zeros(len(logy))
    SSRES=np.zeros(len(logy))

    fx = ajuste( logx, a, b)
    for j in range (0, len(SSTOT)):

        SSTOT[j] = (logy[j] - logymean) ** 2
        SSRES[j] = (logy[j] - fx[j]) ** 2

    sstot = np.sum(SSTOT)
    ssres = np.sum(SSRES)

    Rcuadrado = 1 - (ssres/sstot)

    x = np.zeros(len(logy))
    for i in range (0, len(logy), 1):
        x[i] = (logy[i] - fx[i])**2 / fx[i]
    chi2 = abs(np.sum(x))

    return Rcuadrado, a, b, perr[0], perr[1]

def polinomio (logx, logy):
	n = 3
	def ajuste_p(x,a3, a2, a1, a0):
		return a3 * x**3 + a2 * x**2 + a1 * x + a0

	popt,pcov = curve_fit (ajuste_p, logx, logy)
	perr = np.sqrt(np.diag(pcov))
	a = np.zeros(n + 1)
	da = np.zeros(n + 1)

	for i in range (0, n + 1):
		a[i] = popt [n - i]
		da[i] = perr[n - i]
	logymean = np.mean(logy)
	SSTOT=np.zeros(len(logy))
	SSRES=np.zeros(len(logy))
	fx = ajuste_p(logx, a[3], a[2], a[1], a[0])
	for j in range (0, len(SSTOT)):
		SSTOT[j] = (logy[j] - logymean) ** 2
		SSRES[j] = (logy[j] - fx[j]) ** 2

	sstot = np.sum(SSTOT)
	ssres = np.sum(SSRES)

	Rcuadrado = 1 - (ssres/sstot)

	x = np.zeros(len(logy))
	for i in range (0, len(logy), 1):
		x[i] = (logy[i] - fx[i])**2 / fx[i]
	chi2 = abs(np.sum(x))

	return Rcuadrado, a, da

M = np.loadtxt('barrido_rho_0.13.txt', delimiter = '\t')
T = M[:, 0]
E = M[:, 1]
aceptacion_x = M[:, 2]
aceptacion_p = M[:, 3]

'''
dx  = M[:, 4]
dp = M[:, 5]
Rp, ap, bp, dap, dbp = ajuste_raiz (T, dp)
x = np.linspace(min(T), max(T), 1000)

print(Rp, '\t', ap, '\t', bp)#, '\t', cp)
'''

M2 = np.loadtxt('barrido_rho_0.14.txt', delimiter = '\t')
T2 = M2[:, 0]
E2 = M2[:, 1]
aceptacion_x2 = M2[:, 2]
aceptacion_p2 = M2[:, 3]

M3 = np.loadtxt('barrido_rho_0.16.txt', delimiter = '\t')
T3 = M3[:, 0]
E3 = M3[:, 1]
aceptacion_x3 = M3[:, 2]
aceptacion_p3 = M3[:, 3]

plt.figure()

plt.grid()
plt.ylabel('E (MeV)', fontsize = 13)
plt.xlabel('T (MeV)', fontsize = 13)
plt.title('Curva calórica', fontsize = 15)

plt.plot(T, E, '.', label = 'rho = 0.13')
plt.plot(T2, E2, '.', label = 'rho = 0.14')
plt.plot(T3, E3, '.', label = 'rho = 0.16')

#plt.plot(x, np.sqrt(x) * ap + bp, label = 'ajuste dp')
#plt.plot(T, dx, '.', label = 'rho = 0.10')
plt.legend()
plt.show()

