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

def deltas(N):
	for i in range(0, N):
		j = 16 - i
		M = np.loadtxt('deltas_rho_0.%i.txt' %j, delimiter = '\t')
		T = M[:, 0]
		E = M[:, 1] / 512
		aceptacion_x = M[:, 2]
		aceptacion_p = M[:, 3]
		Dx = M[:, 4]
		Dp= M[:, 5]
		R_x, ax, bx, dax, dbx = ajuste_raiz (T, Dx)
		R_p, ap, bp, dap, dbp = ajuste_raiz (T, Dp)
		x = np.linspace(min(T), max(T), 1000)
		plt.plot(x, np.sqrt(x) * ax + bx, label = 'ajuste dx')
		plt.plot(T, Dx, '.', label = 'rho = 0.%i' %j)

		plt.plot(x, np.sqrt(x) * ap + bp, label = 'ajuste dp')
		plt.plot(T, Dp, '.', label = 'rho = 0.%i' %j)
		print('rho ---> 0.', j)
#		print('R^2_x \t | \t ax \t | \t bx \t | \t R^2_p \t | \t ap \t | \t bp')
#		print(j / 100.0, '\t', R_x, '\t', ax, '\t', bx, '\t', R_p, '\t', ap, '\t', bp)
		print(j / 100.0, '\t double ax = ', ax, ', bx = ', bx, ', ap = ', ap, ', bp = ', bp, ';')
	return T, E, aceptacion_x, aceptacion_p, Dx, Dp
'''
dx  = M[:, 4]
dp = M[:, 5]
Rp, ap, bp, dap, dbp = ajuste_raiz (T, dp)
x = np.linspace(min(T), max(T), 1000)

M2 = np.loadtxt('5barrido_rho_0.16.txt', delimiter = '\t')
T2 = M2[:, 0]
E2 = M2[:, 1]
#aceptacion_x2 = M2[:, 2]
#aceptacion_p2 = M2[:, 3]
print(Rp, '\t', ap, '\t', bp)#, '\t', cp)
'''
plt.figure()

plt.grid()
plt.ylabel('E (MeV)', fontsize = 13)
plt.xlabel('T (MeV)', fontsize = 13)
plt.title('Deltas para Metrópolis', fontsize = 15)

T, E, a_c, a_p, dx, dp = deltas(5)

plt.legend()
plt.show()

