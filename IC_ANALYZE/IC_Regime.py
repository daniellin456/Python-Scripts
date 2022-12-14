import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import *


def cooling(nH, T, x):
    G0 = 1.0 / 1.7
    bet = 0.74 / (T ** 0.068)
    ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
    if x >= 0.1:
        param = G0 * np.sqrt(T) / (nH * 0.1)
        CII_Cooling = 92.0 * 1.38e-16 * 2.0 * (
                (2.8E-7 * (T / 100.0) ** -0.5) * 0.1 + 8e-10 * (T / 100.0) ** 0.07) * 3.5e-4 * 0.4 * np.exp(
            -92.0 / T)
        OI_Cooling = 1e-26 * np.sqrt(T) * (24.0 * np.exp(-228.0 / T) + 7.0 * np.exp(-326.0 / T)) * 4.5E-4
        H2_Cooling = 7.3E-19 * 0.1 * np.exp(-118400.0 / T)
        CII_Line_Cooling = 6.2e4 * 1.38e-16 * 1.0 * (2.3e-8 * (T / 10000.0) ** -0.5 * 0.1 + 1e-12) * np.exp(
            -6.2e4 / T) * 3.5e-4 * 0.4
        OI_Line_Cooling = 4.5e-4 * (2.3e4 * 1.38e-16 / 3.0 * (5.1e-9 * (T / 10000.) ** 0.57 * 0.1 + 1e-12) * np.exp(
            -2.3e4 / T) + 4.9e4 * 1.38e-16 / 3.0 * (2.5e-9 * (T / 10000.) ** 0.57 * 0.1 + 1e-12) * np.exp(
            -4.9e4 / T) + 2.6e4 * 1.38e-16 * 1.0 * (5.2e-9 * (T / 10000.) ** 0.57 * 0.1 + 1e-12) * np.exp(-2.6e4 / T))
        Recombination_Cooling = 4.65E-30 * (T ** 0.94) * (param ** bet) * 0.1
    elif x <= 1.4e-4:
        param = G0 * np.sqrt(T) / (nH * 1.4e-4)
        CII_Cooling = 92.0 * 1.38e-16 * 2.0 * (
                (2.8E-7 * (T / 100.0) ** -0.5) * 1.4e-4 + 8e-10 * (T / 100.0) ** 0.07) * 3.5e-4 * 0.4 * np.exp(
            -92.0 / T)
        OI_Cooling = 1e-26 * np.sqrt(T) * (24.0 * np.exp(-228.0 / T) + 7.0 * np.exp(-326.0 / T)) * 4.5E-4
        H2_Cooling = 7.3E-19 * 1.4e-4 * np.exp(-118400.0 / T)
        CII_Line_Cooling = 6.2e4 * 1.38e-16 * 1.0 * (2.3e-8 * (T / 10000.0) ** -0.5 * 1.4e-4 + 1e-12) * np.exp(
            -6.2e4 / T) * 3.5e-4 * 0.4
        OI_Line_Cooling = 4.5e-4 * (2.3e4 * 1.38e-16 / 3.0 * (5.1e-9 * (T / 10000.) ** 0.57 * 1.4e-4 + 1e-12) * np.exp(
            -2.3e4 / T) + 4.9e4 * 1.38e-16 / 3.0 * (2.5e-9 * (T / 10000.) ** 0.57 * 1.4e-4 + 1e-12) * np.exp(
            -4.9e4 / T) + 2.6e4 * 1.38e-16 * 1.0 * (5.2e-9 * (T / 10000.) ** 0.57 * 1.4e-4 + 1e-12) * np.exp(
            -2.6e4 / T))
        Recombination_Cooling = 4.65E-30 * (T ** 0.94) * (param ** bet) * 1.4e-4

    else:
        param = G0 * np.sqrt(T) / (nH * ne / nH)
        CII_Cooling = 92.0 * 1.38e-16 * 2.0 * (
                (2.8E-7 * (T / 100.0) ** -0.5) * ne / nH + 8e-10 * (T / 100.0) ** 0.07) * 3.5e-4 * 0.4 * np.exp(
            -92.0 / T)
        OI_Cooling = 1e-26 * np.sqrt(T) * (24.0 * np.exp(-228.0 / T) + 7.0 * np.exp(-326.0 / T)) * 4.5E-4
        H2_Cooling = 7.3E-19 * ne / nH * np.exp(-118400.0 / T)
        CII_Line_Cooling = 6.2e4 * 1.38e-16 * 1.0 * (2.3e-8 * (T / 10000.0) ** -0.5 * ne / nH + 1e-12) * np.exp(
            -6.2e4 / T) * 3.5e-4 * 0.4
        OI_Line_Cooling = 4.5e-4 * (2.3e4 * 1.38e-16 / 3.0 * (5.1e-9 * (T / 10000.) ** 0.57 * ne / nH + 1e-12) * np.exp(
            -2.3e4 / T) + 4.9e4 * 1.38e-16 / 3.0 * (2.5e-9 * (T / 10000.) ** 0.57 * ne / nH + 1e-12) * np.exp(
            -4.9e4 / T) + 2.6e4 * 1.38e-16 * 1.0 * (5.2e-9 * (T / 10000.) ** 0.57 * ne / nH + 1e-12) * np.exp(
            -2.6e4 / T))
        Recombination_Cooling = 4.65E-30 * (T ** 0.94) * (param ** bet) * ne / nH
    Total_Cooling = CII_Cooling + OI_Cooling + H2_Cooling + OI_Line_Cooling + CII_Line_Cooling + Recombination_Cooling
    return Total_Cooling


def heating(nH, T, x):
    G0 = 1.0 / 1.7
    ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
    if x >= 0.1:
        param = G0 * np.sqrt(T) / (nH * 0.1)
        epsilon = 4.9E-2 / (1.0 + (param / 1925.) ** 0.73) + 3.7E-2 * (T / 1e4) ** 0.7 / (1. + (param / 5e3))
        UV_Heating = 1e-24 * epsilon * G0
    elif x <= 1.4e-4:
        param = G0 * np.sqrt(T) / (nH * 1.4e-4)
        epsilon = 4.9E-2 / (1.0 + (param / 1925.) ** 0.73) + 3.7E-2 * (T / 1e4) ** 0.7 / (1. + (param / 5e3))
        UV_Heating = 1e-24 * epsilon * G0
    else:
        param = G0 * np.sqrt(T) / (nH * ne / nH)
        epsilon = 4.9E-2 / (1.0 + (param / 1925.) ** 0.73) + 3.7E-2 * (T / 1e4) ** 0.7 / (1. + (param / 5e3))
        UV_Heating = 1e-24 * epsilon * G0
    Total_Heating = UV_Heating
    return Total_Heating


def dCoolingdRho(nH, T, x):
    G0 = 1.0 / 1.7
    bet = 0.74 / (T ** 0.068)
    ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
    if x >= 0.1:
        dCII_CoolingdRho = 0
        dOI_CoolingdRho = 0
        dH2_CoolingdRho = 0
        dCII_Line_CoolingdRho = 0
        dOI_Line_CoolingdRho = 0
        dRecombinationdRho = -3.441e-31 * T ** 0.872 * (5.88235294117647 * T ** 0.5 / nH) ** bet / nH
    elif x <= 1.4e-4:
        dCII_CoolingdRho = 0
        dOI_CoolingdRho = 0
        dH2_CoolingdRho = 0
        dCII_Line_CoolingdRho = 0
        dOI_Line_CoolingdRho = 0
        dRecombinationdRho = -4.8174e-34 * T ** 0.872 * (4201.68067226891 * T ** 0.5 / nH) ** bet / nH
    else:
        dCII_CoolingdRho = -1.51085996659307e-26 * T ** -0.25 * np.exp(-92.0 / T) * nH ** -2
        dOI_CoolingdRho = 0
        dH2_CoolingdRho = -1.108062092123e-21 * T ** 0.25 * np.exp(-118400.0 / T) * nH ** -2
        dCII_Line_CoolingdRho = -4.18184455039153e-24 * T ** -0.25 * np.exp(-6.2e4) * nH ** -2
        dOI_Line_CoolingdRho = -1.93423659159274e-29 * T ** 0.82 * np.exp(
            -2.3e4 / T) * nH ** -2 - 2.0199828002567e-29 * T ** 0.82 * np.exp(
            -4.9e4 / T) * nH ** -2 - 6.68820427578872e-29 * T ** 0.82 * np.exp(-2.6e4 / T) * nH ** -2
        dRecombinationdRho = -7.05820373749582e-33 * T ** 1.19 * (387.534026981419 * T ** 0.25) ** bet

    dTotal_CoolingdRho = dCII_CoolingdRho + dOI_CoolingdRho + dH2_CoolingdRho + dCII_Line_CoolingdRho + dOI_Line_CoolingdRho + dRecombinationdRho
    return dTotal_CoolingdRho


def dHeatingdRho(nH, T, x):
    if x >= 0.1:
        dUV_HeatingdRho = 4.05820402221877e-32 * T ** 1.2 / (
                nH ** 2 * (0.00117647058823529 * T ** 0.5 / nH + 1.0) ** 2) + 3.07049129489545e-28 * (
                                  T ** 0.5 / nH) ** 0.73 / (
                                  nH * (0.0145927738365174 * (T ** 0.5 / nH) ** 0.73 + 1.0) ** 2)
    elif x <= 1.4e-4:
        dUV_HeatingdRho = 2.89871715872769e-29 * T ** 1.2 / (
                nH ** 2 * (0.840336134453782 * T ** 0.5 / nH + 1.0) ** 2) + 1.19016077589066e-26 * (
                                  T ** 0.5 / nH) ** 0.73 / (nH * ((T ** 0.5 / nH) ** 0.73 + 0.565634140065452) ** 2)
    else:
        dUV_HeatingdRho = 0
    return dUV_HeatingdRho


def dCoolingdTemp(nH, T, x):
    G0 = 1.0 / 1.7
    bet = 0.74 / (T ** 0.068)
    ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
    if x >= 0.1:
        param = G0 * np.sqrt(T) / (nH * 0.1)
        dCII_CoolingdTemp = (-4.976832e-25 / T ** 1.5 + 1.44215842722077e-28 / T ** 0.93) * np.exp(-92.0 / T) + 92.0 * (
                9.953664e-25 / T ** 0.5 + 2.0602263246011e-27 * T ** 0.07) * np.exp(-92.0 / T) / T ** 2
        dOI_CoolingdTemp = 2.25e-30 * (
                7.0 * np.exp(-326.0 / T) + 24.0 * np.exp(-228.0 / T)) / T ** 0.5 + 4.5e-30 * T ** 0.5 * (
                                   2282.0 * np.exp(-326.0 / T) / T ** 2 + 5472.0 * np.exp(-228.0 / T) / T ** 2)
        dH2_CoolingdTemp = 8.6432e-15 * np.exp(-118400.0 / T) / T ** 2
        dCII_Line_CoolingdTemp = -1.377516e-22 * np.exp(-62000.0 / T) / T ** 1.5 + 8.68 * (
                1.96788e-18 / T ** 0.5 + 8.556e-24) * np.exp(-62000.0 / T) / T ** 2
        dRecombination_CoolingdTemp = 4.371e-31 * (5.88235294117647 * T ** 0.5 / nH) ** (
                0.74 / T ** 0.068) / T ** 0.0600000000000001 + 4.65e-31 * T ** 0.94 * (
                                              5.88235294117647 * T ** 0.5 / nH) ** (0.74 / T ** 0.068) * (
                                              -0.05032 * np.log(
                                          5.88235294117647 * T ** 0.5 / nH) / T ** 1.068 + 0.37 / T ** 1.068)
    elif x <= 1.4e-4:
        param = G0 * np.sqrt(T) / (nH * 1.4e-4)
        dCII_CoolingdTemp = (-6.9675648e-28 / T ** 1.5 + 1.44215842722077e-28 / T ** 0.93) * np.exp(
            -92.0 / T) + 92.0 * (1.39351296e-27 / T ** 0.5 + 2.0602263246011e-27 * T ** 0.07) * np.exp(
            -92.0 / T) / T ** 2
        dOI_CoolingdTemp = 2.25e-30 * (
                7.0 * np.exp(-326.0 / T) + 24.0 * np.exp(-228.0 / T)) / T ** 0.5 + 4.5e-30 * T ** 0.5 * (
                                   2282.0 * np.exp(-326.0 / T) / T ** 2 + 5472.0 * np.exp(-228.0 / T) / T ** 2)
        dH2_CoolingdTemp = 1.210048e-17 * np.exp(-118400.0 / T) / T ** 2
        dCII_Line_CoolingdTemp = -1.9285224e-25 * np.exp(-62000.0 / T) / T ** 1.5 + 8.68 * (
                2.755032e-21 / T ** 0.5 + 8.556e-24) * np.exp(-62000.0 / T) / T ** 2
        dRecombination_CoolingdTemp = 6.1194e-34 * (4201.68067226891 * T ** 0.5 / nH) ** (
                0.74 / T ** 0.068) / T ** 0.0600000000000001 + 6.51e-34 * T ** 0.94 * (
                                              4201.68067226891 * T ** 0.5 / nH) ** (0.74 / T ** 0.068) * (
                                              -0.05032 * np.log(
                                          4201.68067226891 * T ** 0.5 / nH) / T ** 1.068 + 0.37 / T ** 1.068)
    else:
        param = G0 * np.sqrt(T) / (nH * ne / nH)
        dCII_CoolingdTemp = (-3.77714991648267e-27 / (T ** 1.25 * nH) + 1.44215842722077e-28 / T ** 0.93) * np.exp(
            -92.0 / T) + 92.0 * (1.51085996659307e-26 / (T ** 0.25 * nH) + 2.0602263246011e-27 * T ** 0.07) * np.exp(
            -92.0 / T) / T ** 2
        dOI_CoolingdTemp = 2.25e-30 * (
                7.0 * np.exp(-326.0 / T) + 24.0 * np.exp(-228.0 / T)) / T ** 0.5 + 4.5e-30 * T ** 0.5 * (
                                   2282.0 * np.exp(-326.0 / T) / T ** 2 + 5472.0 * np.exp(-228.0 / T) / T ** 2)
        dH2_CoolingdTemp = 1.31194551707363e-16 * np.exp(-118400.0 / T) / (
                T ** 1.75 * nH) + 2.7701552303075e-22 * np.exp(-118400.0 / T) / (T ** 0.75 * nH)
        dCII_Line_CoolingdTemp = -1.04546113759788e-24 * np.exp(-62000.0 / T) / (T ** 1.25 * nH) + 8.68 * (
                2.98703182170823e-20 / (T ** 0.25 * nH) + 8.556e-24) * np.exp(-62000.0 / T) / T ** 2
        dRecombination_CoolingdTemp = 8.39926244762003e-33 * T ** 0.19 * (387.534026981419 * T ** 0.25) ** (
                0.74 / T ** 0.068) / nH + 7.05820373749582e-33 * T ** 1.19 * (387.534026981419 * T ** 0.25) ** (
                                              0.74 / T ** 0.068) * (-0.05032 * np.log(
            387.534026981419 * T ** 0.25) / T ** 1.068 + 0.185 / T ** 1.068) / nH

    dTotal_CoolingdTemp = dCII_CoolingdTemp + dOI_CoolingdTemp + dH2_CoolingdTemp + dCII_Line_CoolingdTemp + dRecombination_CoolingdTemp
    return dTotal_CoolingdTemp


def dHeatingdTemp(nH, T, x):
    G0 = 1.0 / 1.7
    ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
    if x >= 0.1:
        param = G0 * np.sqrt(T) / (nH * 0.1)
        dUV_HeatingdTemp = -1.53524564744772e-28 * (T ** 0.5 / nH) ** 0.73 / (
                T ** 1.0 * (0.0145927738365174 * (T ** 0.5 / nH) ** 0.73 + 1.0) ** 2) + 2.41463139322017e-29 / (
                                   T ** 0.3 * (
                                   0.00117647058823529 * T ** 0.5 / nH + 1.0)) - 2.02910201110938e-32 * T ** 0.2 / (
                                   nH * (0.00117647058823529 * T ** 0.5 / nH + 1.0) ** 2)
    elif x <= 1.4e-4:
        param = G0 * np.sqrt(T) / (nH * 1.4e-4)
        dUV_HeatingdTemp = -5.9508038794533e-27 * (T ** 0.5 / nH) ** 0.73 / (
                T ** 1.0 * ((T ** 0.5 / nH) ** 0.73 + 0.565634140065452) ** 2) + 2.41463139322017e-29 / (
                                   T ** 0.3 * (
                                   0.840336134453782 * T ** 0.5 / nH + 1.0)) - 1.44935857936385e-29 * T ** 0.2 / (
                                   nH * (0.840336134453782 * T ** 0.5 / nH + 1.0) ** 2)
    else:
        param = G0 * np.sqrt(T) / (nH * ne / nH)
        dUV_HeatingdTemp = -1.63245709529796e-27 / (
                T ** 0.8175 * (0.310335707241436 * T ** 0.1825 + 1.0) ** 2) + 2.41463139322017e-29 / (
                                   T ** 0.3 * (0.0775068053962838 * T ** 0.25 + 1.0)) - 6.68394162493118e-31 / (
                                   T ** 0.05 * (0.0775068053962838 * T ** 0.25 + 1.0) ** 2)
    return dUV_HeatingdTemp


def plot_equilibrium_temperature(nH_array, balance_temperature):
    plt.figure(figsize=(5, 4))
    plt.loglog(nH_array, balance_temperature, linewidth=2, linestyle="-", color="r")
    plt.title("Balance temperature")
    plt.tick_params(which="major", width=1, length=5)
    plt.tick_params(which="minor", width=1, length=3)
    plt.plot()
    plt.show()
    return


def calculate_equilibrium_temperature(nH_array, T_array):
    total_cooling_rate = np.zeros(shape=len(T_array))
    balance_temperature = np.zeros(shape=len(nH_array))

    for i in range(0, len(nH_array)):
        nH = nH_array[i]

        for j in range(0, len(T_array)):
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH

            total_cooling_rate[j] = np.abs(nH ** 2 * cooling(nH, T, x) - nH * heating(nH, T, x))

        index = np.argsort(total_cooling_rate)[0]
        balance_temperature[i] = T_array[index]

    plot_equilibrium_temperature(nH_array, balance_temperature)
    return


def calculate_Gamma(m, n):
    gamma = 5 / 3
    c = 2 * (3 - 2 * m) / (5 - 2 * m - 2 * n)
    d = ((3 - 2 * m) / (5 - 2 * m - 2 * n) - 1) * 2
    Gamma1 = (gamma * 2 * (d / 3 + 1) + (n - m) * (c - (gamma - 1) * d)) / (
            2 * (d / 3 + 1) + (n - 1) * (c - (gamma - 1) * d))
    Gamma2 = 1 + (3 - 2 * m) / (2 * (n - 1))
    Gamma3 = 1 + (2 - m) / (n - 1)
    return Gamma1, Gamma2, Gamma3


def plot_m(fig, ax, extent, m_matrix):
    print("Extrema of m: min:%24.14e, max: %24.14e" % (np.min(m_matrix), np.max(m_matrix)))

    im = ax.imshow(m_matrix.T, origin='lower', cmap="rainbow", norm=LogNorm(vmin=1e-1, vmax=1e2), extent=extent)
    ax.set_title("m(nH, T)")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend='both', label='m')
    return


def plot_n(fig, ax, extent, n_matrix):
    print("Extrema of n: min:%24.14e, max: %24.14e" % (np.min(n_matrix), np.max(n_matrix)))

    im = ax.imshow(n_matrix.T, origin='lower', cmap="rainbow", norm=LogNorm(vmin=1e-2, vmax=1e2), extent=extent)
    ax.set_title("n(nH, T)")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend='both', label='n')
    return


def plot_Gamma1(fig, ax, extent, Gamma1):
    print("Extrema of Gamma1: min:%24.14e, max: %24.14e" % (np.min(Gamma1), np.max(Gamma1)))

    im = ax.imshow(Gamma1.T, origin='lower', cmap="rainbow", norm=SymLogNorm(linthresh=1, vmin=-10, vmax=10),
                   extent=extent)
    ax.set_title(r"$\Gamma_1\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_1$")
    return


def plot_Gamma2(fig, ax, extent, Gamma2):
    print("Extrema of Gamma2: min:%24.14e, max: %24.14e" % (np.min(Gamma2), np.max(Gamma2)))

    im = ax.imshow(Gamma2.T, origin="lower", cmap="rainbow", norm=SymLogNorm(linthresh=1, vmin=-10, vmax=10),
                   extent=extent)
    ax.set_title(r"$\Gamma_2\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_2$")
    return


def plot_Gamma3(fig, ax, extent, Gamma3):
    print("Extrema of Gamma_3: min: %24.14e, max: %24.14e" % (np.min(Gamma3), np.max(Gamma3)))

    im = ax.imshow(Gamma3.T, origin="lower", cmap="rainbow", norm=SymLogNorm(linthresh=0.01, vmin=-10, vmax=10),
                   extent=extent)
    ax.set_title(r"$\Gamma_3\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_3$")
    return


def plot_Gamma1_m1(fig, ax, extent, Gamma1_m1):
    return


def plot_Gamma2_m1(fig, ax, extent, Gamma2_m1):
    return


def plot_Gamma3_m1(fig, ax, extent, Gamma3_m1):
    return


def make_plot(extent, m_matrix, n_matrix, Gamma1, Gamma2, Gamma3):
    fig, ((ax0, ax1, ax2), (ax3, ax4, ax5)) = plt.subplots(nrows=2, ncols=3, figsize=(21, 12))

    plot_m(fig, ax0, extent, m_matrix)
    plot_n(fig, ax1, extent, n_matrix)
    plot_Gamma1(fig, ax3, extent, Gamma1)
    plot_Gamma2(fig, ax4, extent, Gamma2)
    plot_Gamma3(fig, ax5, extent, Gamma3)

    ## quiver plot
    # plot_Gamma1_m1(fig, ax3, nH_list, T_list, Gamma1 - 1)
    # plot_Gamma2_m1(fig, ax4, nH_list, T_list, Gamma2 - 1)
    # plot_Gamma3_m1(fig, ax5, nH_list, T_list, Gamma3 - 1)

    plt.tight_layout(pad=2.5)
    plt.show()
    # plt.savfig("IC_Analyze.png", dpi=600)
    return


def mesh_method():
    nH_start = -2
    nH_end = 6
    nH_count = 100
    T_start = -2
    T_end = 6
    T_count = 100

    nH_array = np.logspace(nH_start, nH_end, nH_count)
    T_array = np.logspace(T_start, T_end, T_count)
    m_matrix = np.zeros(shape=(nH_count, T_count))
    n_matrix = np.zeros(shape=(nH_count, T_count))

    nH_mesh, T_mesh = np.meshgrid(nH_array, T_array);

    for i in range(0, nH_count):
        for j in range(0, T_count):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            total_cooling_rate = np.abs(nH ** 2 * cooling(nH, T, x) - nH * heating(nH, T, x))
            m_matrix[i, j] = nH / total_cooling_rate * np.abs(
                heating(nH, T, x) + nH * dHeatingdRho(nH, T, x) - 2 * nH * cooling(nH, T, x) - nH ** 2 * dCoolingdRho(
                    nH, T, x))
            n_matrix[i, j] = T / total_cooling_rate * np.abs(
                nH * dHeatingdTemp(nH, T, x) - nH ** 2 * dCoolingdTemp(nH, T, x))

    plt.title("m(nH, T)")
    sc = plt.scatter(nH_mesh, T_mesh, c=m_matrix.T, cmap="rainbow", norm=LogNorm(vmin=1e-1, vmax=1e2))
    plt.xscale("log")
    plt.yscale("log")
    plt.colorbar(sc)
    plt.show()

    plt.title("n(nH, T)")
    sc = plt.scatter(nH_mesh, T_mesh, c=n_matrix.T, cmap="rainbow", norm=LogNorm(vmin=1e-2, vmax=1e2))
    plt.xscale("log")
    plt.yscale("log")
    plt.colorbar(sc)
    plt.show()
    # Gamma1, Gamma2, Gamma3 = calculate_Gamma(m_matrix, n_matrix)
    # make_plot(nH_array, T_array, m_matrix, n_matrix, Gamma1, Gamma2, Gamma3)


def main():
    nH_start = -2
    nH_end = 6
    nH_count = 100
    T_start = -2
    T_end = 6
    T_count = 100
    extent = [nH_start, nH_end, T_start, T_end]

    nH_array = np.logspace(nH_start, nH_end, nH_count)
    T_array = np.logspace(T_start, T_end, T_count)
    m_matrix = np.zeros(shape=(nH_count, T_count))
    n_matrix = np.zeros(shape=(nH_count, T_count))
    # calculate_equilibrium_temperature(nH_array, T_array)

    for i in range(0, nH_count):
        for j in range(0, T_count):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH
            total_cooling_rate = np.abs(nH ** 2 * cooling(nH, T, x) - nH * heating(nH, T, x))
            m_matrix[i, j] = nH / total_cooling_rate * np.abs(
                heating(nH, T, x) + nH * dHeatingdRho(nH, T, x) - 2 * nH * cooling(nH, T, x) - nH ** 2 * dCoolingdRho(
                    nH, T, x))
            n_matrix[i, j] = T / total_cooling_rate * np.abs(
                nH * dHeatingdTemp(nH, T, x) - nH ** 2 * dCoolingdTemp(nH, T, x))

    Gamma1, Gamma2, Gamma3 = calculate_Gamma(m_matrix, n_matrix)
    make_plot(extent, m_matrix, n_matrix, Gamma1, Gamma2, Gamma3)


if __name__ == "__main__":
    main()
    mesh_method();
