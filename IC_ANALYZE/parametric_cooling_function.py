import numpy as np

kb = 1.3606e-16
gamma = 5 / 3
G = 6.67e-8
mu = 2.31
mp = 1.66e-24


def main():
    nH = 1280  # density (cm^-3)
    T = 80  # temperature (k)
    m = 1.3  # power index of density
    n = 2.3  # power index of temperature
    alpha = 1  # ratio of cooling time to free fall time

    lambda_0 = kb / (gamma - 1) * np.sqrt(32 * G * mu * mp / 3 / np.pi) * nH ** (1.5 - m) * T ** (1 - n) / alpha
    free_fall_time = np.sqrt(3 * np.pi / 32 / G / mu / mp / nH)
    cooling_time = nH * kb * T / (gamma - 1) / nH ** m / T ** n

    lambda_1 = cooling_time / free_fall_time / alpha
    print("  density (cm^-3) =", nH)
    print("  temperature (k) =", T)
    print("  m               =", m)
    print("  n               =", n)
    print("  t_cooling/t_ff  =", alpha)
    print("  lambda_0        =", lambda_0)
    print("  lambda_1        =", lambda_1)

    nH = 128000
    T = 80
    m = 1.3
    n = 2.3
    cooling_time = nH * kb * T / (gamma - 1) / nH ** m / T ** n / lambda_0
    free_fall_time = np.sqrt(3 * np.pi / 32 / G / mu / mp / nH)
    print("  Ratio           =", cooling_time / free_fall_time)

    nH = 1280
    T = 400
    m = 1.3
    n = 2.3
    cooling_time = nH * kb * T / (gamma - 1) / nH ** m / T ** n / lambda_0
    free_fall_time = np.sqrt(3 * np.pi / 32 / G / mu / mp / nH)
    print("  Ratio           =", cooling_time / free_fall_time)

    return


if __name__ == "__main__":
    main()
