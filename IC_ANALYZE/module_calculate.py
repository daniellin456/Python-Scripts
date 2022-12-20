from Cooling_Function import *


def calculate_total_cooling_rate(nH, T, x):
    return nH * heating(nH, T, x) - nH ** 2 * cooling(nH, T, x)


def calculate_equilibrium_temperature(nH_array, T_array):
    total_cooling_rate = np.zeros(shape=len(T_array))
    balance_temperature = np.zeros(shape=len(nH_array))

    for i in range(0, len(nH_array)):
        nH = nH_array[i]

        for j in range(0, len(T_array)):
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH

            total_cooling_rate[j] = np.abs(calculate_total_cooling_rate(nH, T, x))

        index = np.argsort(total_cooling_rate)[0]
        balance_temperature[i] = T_array[index]

    return balance_temperature


def calculate_m(nH_array, T_array, nH_count, T_count):
    m_matrix = np.zeros(shape=(nH_count, T_count))
    for i in range(0, nH_count):
        for j in range(0, T_count):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            total_cooling_rate = np.abs(calculate_total_cooling_rate(nH, T, x))
            m_matrix[i, j] = nH / total_cooling_rate * np.abs(
                heating(nH, T, x) + nH * dHeatingdRho(nH, T, x) - 2 * nH * cooling(nH, T, x) - nH ** 2 * dCoolingdRho(
                    nH, T, x))
    return m_matrix


def calculate_n(nH_array, T_array, nH_count, T_count):
    n_matrix = np.zeros(shape=(nH_count, T_count))
    for i in range(0, nH_count):
        for j in range(0, T_count):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            total_cooling_rate = np.abs(calculate_total_cooling_rate(nH, T, x))

            n_matrix[i, j] = T / total_cooling_rate * np.abs(
                nH * dHeatingdTemp(nH, T, x) - nH ** 2 * dCoolingdTemp(nH, T, x))
    return n_matrix


def calculate_Gamma(m, n):
    gamma = 5 / 3
    c = 2 * (3 - 2 * m) / (5 - 2 * m - 2 * n)
    d = ((3 - 2 * m) / (5 - 2 * m - 2 * n) - 1) * 2
    Gamma1 = (gamma * 2 * (d / 3 + 1) + (n - m) * (c - (gamma - 1) * d)) / (
            2 * (d / 3 + 1) + (n - 1) * (c - (gamma - 1) * d))
    Gamma2 = 1 + (3 - 2 * m) / (2 * (n - 1))
    Gamma3 = 1 + (2 - m) / (n - 1)
    return Gamma1, Gamma2, Gamma3


def calculate_dLambda_dRho_by_np_diff(nH_array, T_array):
    Lambda = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        nH = nH_array[i]
        for j in range(0, len(T_array)):
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            Lambda[i, j] = calculate_total_cooling_rate(nH, T, x)

    dLambda_dRho = np.zeros(shape=(len(T_array), len(nH_array) - 1))
    for i in range(0, len(T_array)):
        dLambda_dRho[i] = nH_array[1:] / np.abs(Lambda.T[i, 1:]) * np.abs(np.diff(Lambda.T[i]) / np.diff(nH_array))
    return dLambda_dRho


def calculate_dLambda_dRho_by_analytical(nH_array, T_array):
    dLambda_dRho_by_analytical = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            Lambda = calculate_total_cooling_rate(nH, T, x)
            dLambda_dRho_by_analytical[i, j] = nH / np.abs(Lambda) * np.abs(
                heating(nH, T, x) + nH * dHeatingdRho(nH, T, x) - 2 * nH * cooling(nH, T, x) - nH ** 2 * dCoolingdRho(
                    nH, T, x))

    return dLambda_dRho_by_analytical.T


def calculate_dLambda_dTemp_by_np_diff(nH_array, T_array):
    specific_Rho_Array = np.asarray([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e6, 1e6])
    Lambda = np.zeros(shape=(len(specific_Rho_Array), len(T_array)))
    for i in range(0, len(specific_Rho_Array)):
        nH = specific_Rho_Array[i]
        for j in range(0, len(T_array)):
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            Lambda[i, j] = calculate_total_cooling_rate(nH, T, x)

    dLambda_dTemp = np.zeros(shape=(len(specific_Rho_Array), len(T_array) - 1))
    for i in range(0, len(specific_Rho_Array)):
        dLambda_dTemp[i] = np.abs(np.diff(Lambda[i])) / np.diff(T_array) * T_array[1:] / np.abs(
            Lambda[i, 1:])
    return dLambda_dTemp
