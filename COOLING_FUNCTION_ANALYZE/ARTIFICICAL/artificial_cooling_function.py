
from plot_style import *
import numpy as np
import matplotlib.pyplot as plt
plot_style()

m = 1.3
n = 2.3
rho0 = 2 * 1.67 * 1.66e-24 * 1e6
T0 = 10
Lambda0 = 5e-20 / rho0 ** m / T0 ** n

rho = np.logspace( -24, -14, 100 )
T = np.logspace( 1, 3, 100)

Lambda = []
for Temp in T:
	Lambda.append( Lambda0 * rho ** m * Temp ** n )

plt.figure( figsize = (20, 10) )
for i in range( 0, 10, 1 ):
	plt.loglog( rho, Lambda[i * 10], label="T = {:.2f}".format(T[i*10]) )

plt.xlabel(r"$log_{10} \; H_2 \; density (g\;cm^{-3})$")
plt.ylabel(r"$log_{10} \; \Lambda \; (erg\;cm^{-3}\;s^{-1})$")
plt.legend()
plt.plot()
plt.savefig('Artifical_Cooling_Function.pdf', dpi=600)
