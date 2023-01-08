from module_cooling import *


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
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4

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
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
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
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            total_cooling_rate = np.abs(calculate_total_cooling_rate(nH, T, x))

            n_matrix[i, j] = T / total_cooling_rate * np.abs(
                nH * dHeatingdTemp(nH, T, x) - nH ** 2 * dCoolingdTemp(nH, T, x))
    return n_matrix


def calculate_a(m, n):
    a = 1 - (3 - 2 * m) / (5 - 2 * n - 2 * m)
    return a


def calculate_b(m, n):
    b = (3 - 2 * m) / (5 - 2 * n - 2 * m)
    return b


def calculate_c(m, n):
    c = 2 * (3 - 2 * m) / (5 - 2 * n - 2 * m)
    return c


def calculate_d(m, n):
    d = ((3 - 2 * m) / (5 - 2 * m - 2 * n) - 1) * 2
    return d


def calculate_Gamma(m, n):
    gamma = 5 / 3
    # c = 2 * (3 - 2 * m) / (5 - 2 * m - 2 * n)
    # d = ((3 - 2 * m) / (5 - 2 * m - 2 * n) - 1) * 2
    c = calculate_c(m, n)
    d = calculate_d(m, n)
    Gamma1 = (gamma * 2 * (d / 3 + 1) + (n - m) * (c - (gamma - 1) * d)) / (
            2 * (d / 3 + 1) + (n - 1) * (c - (gamma - 1) * d))
    Gamma2 = 1 + (3 - 2 * m) / (2 * (n - 1))
    Gamma3 = 1 + (2 - m) / (n - 1)
    return Gamma1, Gamma2, Gamma3


def calculate_dTotal_Lambda_dRho_by_np_diff(nH_array, T_array):
    Total_Lambda = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        nH = nH_array[i]
        for j in range(0, len(T_array)):
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            Total_Lambda[i, j] = calculate_total_cooling_rate(nH, T, x)

    dTotal_Lambda_dRho_by_np_diff = np.zeros(shape=(len(T_array), len(nH_array) - 1))
    for i in range(0, len(T_array)):
        dTotal_Lambda_dRho_by_np_diff[i] = nH_array[1:] / np.abs(Total_Lambda.T[i, 1:]) * np.abs(
            np.diff(Total_Lambda.T[i]) / np.diff(nH_array))
    return dTotal_Lambda_dRho_by_np_diff


def calculate_dTotal_Lambda_dRho_by_analytical(nH_array, T_array):
    dTotal_Lambda_dRho_by_analytical = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            Lambda = calculate_total_cooling_rate(nH, T, x)
            dTotal_Lambda_dRho_by_analytical[i, j] = nH / np.abs(Lambda) * np.abs(
                heating(nH, T, x) + nH * dHeatingdRho(nH, T, x) - 2 * nH * cooling(nH, T, x) - nH ** 2 * dCoolingdRho(
                    nH, T, x))

    return dTotal_Lambda_dRho_by_analytical.T


def calculate_ionization(nH_array, T_array):
    x = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        nH = nH_array[i]
        for j in range(0, len(T_array)):
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5

            if ne / nH >= 0.1:
                x[i, j] = 0.1
            elif ne / nH <= 1.4e-4:
                x[i, j] = 1.4e-4
            else:
                x[i, j] = ne / nH

    return x.T


def calculate_dGamma_dRho_by_np_diff(nH_array, T_array):
    Gamma = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        nH = nH_array[i]
        for j in range(0, len(T_array)):
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            Gamma[i, j] = heating(nH, T, x)

    dGamma_dRho_by_np_diff = np.zeros(shape=(len(T_array), len(nH_array) - 1))
    for i in range(0, len(T_array)):
        dGamma_dRho_by_np_diff[i] = np.diff(Gamma.T[i]) / np.diff(nH_array)
    return dGamma_dRho_by_np_diff


def calculate_dGamma_dRho_by_analytical(nH_array, T_array):
    dGamma_dRho_by_analytical = np.zeros(shape=(len(nH_array), len(T_array)))

    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            dGamma_dRho_by_analytical[i, j] = dHeatingdRho(nH, T, x)

    return dGamma_dRho_by_analytical.T


def calculate_dLambda_dRho_by_np_diff(nH_array, T_array):
    Lambda = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            Lambda[i, j] = cooling(nH, T, x)

    dLambda_dRho_by_np_diff = np.zeros(shape=(len(T_array), len(nH_array) - 1))
    for i in range(0, len(T_array)):
        dLambda_dRho_by_np_diff[i] = np.diff(Lambda.T[i]) / np.diff(nH_array)

    return dLambda_dRho_by_np_diff


def calculate_dLambda_dRho_by_analytical(nH_array, T_array):
    dLambda_dRho_by_analytical = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            dLambda_dRho_by_analytical[i, j] = dCoolingdRho(nH, T, x)
    return dLambda_dRho_by_analytical.T


def calculate_CII_cooling_rate(nH_array, T_array):
    CII_cooling_rate = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            CII_cooling_rate[i, j] = nH ** 2 * CII_Cooling(nH, T, x)
    return CII_cooling_rate


def calculate_OI_cooling_rate(nH_array, T_array):
    OI_cooling_rate = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            OI_cooling_rate[i, j] = nH ** 2 * OI_Cooling(nH, T, x)
    return OI_cooling_rate


def calculate_H2_cooling_rate(nH_array, T_array):
    H2_cooling_rate = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            H2_cooling_rate[i, j] = nH ** 2 * H2_Cooling(nH, T, x)
    return H2_cooling_rate


def calculate_CII_line_cooling_rate(nH_array, T_array):
    CII_line_cooling_rate = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            CII_line_cooling_rate[i, j] = nH ** 2 * CII_Line_Cooling(nH, T, x)
    return CII_line_cooling_rate


def calculate_OI_line_cooling_rate(nH_array, T_array):
    OI_line_cooling_rate = np.zeros(shape=(len(nH_array), len(T_array)))
    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            OI_line_cooling_rate[i, j] = nH ** 2 * OI_Line_Cooling(nH, T, x)
    return OI_line_cooling_rate


def calculate_Recombination_cooling_rate(nH_array, T_array):
    Recombination_cooling_rate = np.zeros(shape=(len(nH_array), len(T_array)))
    G0 = 1.0 / 1.7
    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4

            bet = 0.74 / (T ** 0.068)
            param = G0 * np.sqrt(T) / (nH * x)
            Recombination_cooling_rate[i, j] = nH ** 2 * Recombination_Cooling(nH, T, x, param, bet)
    return Recombination_cooling_rate


def calculate_UV_heating_rate(nH_array, T_array):
    UV_heating_rate = np.zeros(shape=(len(nH_array), len(T_array)))
    G0 = 1.0 / 1.7
    for i in range(0, len(nH_array)):
        for j in range(0, len(T_array)):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            if x >= 0.1:
                x = 0.1
            elif x <= 1.4e-4:
                x = 1.4e-4
            param = G0 * np.sqrt(T) / (nH * x)
            UV_heating_rate[i, j] = nH * UV_Heating(nH, T, x, param, G0)
    return UV_heating_rate
