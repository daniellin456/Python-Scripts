import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

G = 6.67e-8
mu = 2.31
mp = 1.67e-24
gamma = 5/3
kb = 1.3806e-16
lambda_0 = 2.671e-33
m = 1.3
n = 2.3

def func(T, nH):
    return np.sqrt(3*np.pi / 32.0 / G/ mu / mp / nH) - nH * kb * T / (gamma - 1) / lambda_0 / nH ** m / T ** n

def main():
    nH_start = -2
    nH_end = 10
    nH_count = 1000
    root = []

    for nH in nH_array:
        root.append(fsolve(func, 1.0, nH)[0])

    plt.plot(np.log10(nH_array), np.log10(root))
    plt.show()
    return


if __name__ == "__main__":
    main()
