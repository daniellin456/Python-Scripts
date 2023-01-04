import numpy as np

kb = 1.3606e-16
gamma = 5 / 3
G = 6.67e-8
mu = 2.31
mp = 1.66e-24


def main():
    nH = 1280  # density (cm^-3)
    T = 10  # temperature (k)
    m = 2  # power index of density
    n = 4  # power index of temperature
    alpha = 10  # ratio of cooling time to free fall time
    lambda_0 = kb / (gamma - 1) * np.sqrt(32 * G * mu * mp / 3 / np.pi) * nH ** (1.5 - m) * T ** (1 - n) / alpha
    print("  density (cm^-3) =", nH)
    print("  temperature (k) =", T)
    print("  m               =", m)
    print("  n               =", n)
    print("  t_cooling/t_ff  =", alpha)
    print("  lambda_0        =", lambda_0)
    return


if __name__ == "__main__":
    main()
