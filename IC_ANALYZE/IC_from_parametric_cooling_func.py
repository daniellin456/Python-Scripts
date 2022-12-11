import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

'''
Function	: cooling
Description	: Calculate cooling rate of major coolant in interstellar medium
Note		: 1. x = ne/nH, ne is electron density
			  2. All of cooling functions are from ramses_ism/hydro/cooling_module_frig.f90
			  3. Calculate cooling rate in three part (a) 1.4e-4 < x < 0.1, (b) x <= 1.4e-4 (c) x >= 0.1
Parameter	: 1. nH: Number density
			  2. T: Temperature
			  3. x: Analytic ionisation in absence of photoionisation.
Return		: Total cooling rate
'''


def cooling(nH, T, x):
    G0 = 1 / 1.7
    power = 0.74 / (T ** 0.068)

    if x >= 0.1:
        CII = 9.954e-25 * T ** -0.5 * np.exp(-92 / T) + 2.06e-27 * T ** 0.07 * np.exp(-92 / T)
        OI = 1.08e-28 * T ** 0.5 * np.exp(-228 / T) + 3.15e-29 * T ** 0.5 * np.exp(-326 / T)
        H2 = 7.3e-20 * np.exp(-118400 / T)
        CII_M = 2.755e-22 * T ** -0.5 * np.exp(-6.2e4 / T) + 1.198e-27 * np.exp(-6.2e4 / T)
        OI_M = 1.274e-27 * T ** 0.57 * np.exp(-2.3e4 / T) + 4.761e-28 * np.exp(-2.3e4 / T) \
               + 1.331e-27 * T ** 0.57 * np.exp(-4.9e4 / T) + 1.014e-27 * np.exp(-4.9e4 / T) \
               + 4.406e-27 * T ** 0.57 * np.exp(-2.6e4 / T) + 1.615e-27 * np.exp(-2.6e4 / T)
        Rec = 4.65e-31 * T ** 0.94 * (G0 / 0.1 * T ** 0.5 / nH) ** power
    elif x <= 1.4e-4:
        CII = 1.394e-27 * T ** -0.5 * np.exp(-92 / T) + 2.06e-27 * T ** 0.07 * np.exp(-92 / T)
        OI = 1.08e-28 * T ** 0.5 * np.exp(-228 / T) + 3.15e-29 * T ** 0.5 * np.exp(-326 / T)
        H2 = 1.022e-22 * np.exp(-118400 / T)
        CII_M = 3.857e-25 * T ** -0.5 * np.exp(-6.2e4 / T) + 1.198e-27 * np.exp(-6.2e4 / T)
        OI_M = 1.789e-30 * T ** 0.57 * np.exp(-2.3e4 / T) + 4.761e-28 * np.exp(-2.3e4 / T) \
               + 1.863e-30 * T ** 0.57 * np.exp(-4.9e4 / T) + 1.014e-27 * np.exp(-4.9e4 / T) \
               + 6.165e-30 * T ** 0.57 * np.exp(-2.6e4 / T) + 1.615e-27 * np.exp(-2.6e4 / T)
        Rec = 6.51e-34 * T ** 0.94 * (G0 / 1.4e-4 * T ** 0.5 / nH) ** power
    else:
        CII = 1.511e-26 * T ** -0.25 * np.exp(-92 / T) + 2.06e-27 * T ** 0.07 * np.exp(-92 / T)
        OI = 1.08e-28 * T ** 0.5 * np.exp(-228 / T) + 3.15e-29 * T ** 0.5 * np.exp(-326 / T)
        H2 = 1.108e-21 * T ** 0.25 * np.exp(-118400 / T) / nH
        CII_M = 4.182e-24 * T ** -0.25 * np.exp(-6.2e4 / T) / nH + 1.198e-27 * np.exp(-6.2e4 / T)
        OI_M = 1.935e-29 * T ** 0.82 * np.exp(-2.3e4 / T) / nH + 4.761e-28 * np.exp(-2.3e4 / T) \
               + 2.020e-29 * T ** 0.82 * np.exp(-4.9e4 / T) / nH + 1.014e-27 * np.exp(-4.9e4 / T) \
               + 6.687e-29 * T ** 0.82 * np.exp(-2.6e4 / T) / nH + 1.615e-27 * np.exp(-2.6e4 / T)
        Rec = 7.059e-33 * T ** 1.19 * (G0 / 1.518e-3 * T ** 0.25) ** power / nH

    return CII + OI + H2 + CII_M + OI_M + Rec


'''
Function	: heating
Description	: Calculate heating rate of background UV radiation
Note		: 1. x = ne/nH, ne is electron density
			  2. All of heating functions are from ramses_ism/hydro/cooling_module_frig.f90
			  3. Calculate heating rate in three part (a) 1.4e-4 < x < 0.1, (b) x <= 1.4e-4 (c) x >= 0.1
Parameter	: 1. nH: Number density
			  2. T: Temperature
			  3. x: Analytic ionisation in absence of photoionisation.
Return		: Total heating rate
'''


def heating(nH, T, x):
    G0 = 1 / 1.7

    if x >= 0.1:
        UV = G0 * (4.9e-26 / (1 + G0 ** 0.73 * T ** 0.365 * nH ** -0.73 * 0.021)
                   + 5.864e-29 * T ** 0.7 / (1 + G0 * T ** 0.5 * 0.002 / nH))
    elif x <= 1.4e-4:
        UV = G0 * (4.9e-26 / (1 + G0 ** 0.73 * T ** 0.365 * nH ** -0.73 * 2.604)
                   + 5.864e-29 * T ** 0.7 / (1 + G0 * T ** 0.5 * 1.42 / nH))
    else:
        UV = G0 * (4.9e-26 / (1 + G0 ** 0.73 * T ** 0.1825 * nH ** 0.73 * 0.457)
                   + 5.864e-29 * T ** 0.7 / (1 + G0 * T ** 0.5 * 0.132))
    return UV


def dcoolingdrho(nH, T, x):
    G0 = 1 / 1.7
    power = 0.74 / (T ** 0.068)
    if x >= 0.1:
        dCIIdrho = 0
        dOIdrho = 0
        dH2drho = 0
        dCII_mdrho = 0
        dOI_mdrho = 0
        dRecdrho = -3.441e-31 * T ** 0.872 * (G0 / 0.1 * T ** 0.5) ** power * nH ** (power - 1)
    elif x <= 1.4e-4:
        dCIIdrho = 0
        dCIIdrho = 0
        dOIdrho = 0
        dH2drho = 0
        dCII_mdrho = 0
        dOI_mdrho = 0
        dRecdrho = -4.8174e-34 * T ** 0.872 * (G0 / 1.4e-4 * T ** 0.5) ** power * nH ** (power - 1)
    else:
        dCIIdrho = -1.511e-26 * T ** -0.25 * np.exp(-92 / T) * nH ** -2
        dOIdrho = 0
        dH2drho = -1.108e-21 * T ** 0.25 * np.exp(-118400 / T) * nH ** -2
        dCII_mdrho = -4.182e-24 * T ** -0.25 * np.exp(-6.2e4 / T) * nH ** -2
        dOI_mdrho = -1.935e-29 * T ** 0.82 * np.exp(-2.3e4 / T) * nH ** -2 - \
                    2.020e-29 * T ** 0.82 * np.exp(-4.9e4 / T) * nH ** -2 - \
                    6.688e-29 * T ** 0.82 * np.exp(-2.6e4 / T) * nH ** -2
        dRecdrho = -7.059e-33 * T ** 1.19 * nH ** -2 * (G0 / 1.518e-3 * T ** 0.25) ** power

    return dCIIdrho + dOIdrho + dH2drho + dCII_mdrho + dOI_mdrho + dRecdrho


def dcoolingdT(nH, T, x):
    G0 = 1 / 1.7
    power = 0.74 / (T ** 0.068)
    if x >= 0.1:
        dCIIdT = 9.954e-25 * np.exp(-92 / T) * (-0.5 * T ** -1.5 + 92 * T ** -2.5)
        dOIdT = 1.08e-28 * np.exp(-228 / T) * (0.5 * T ** -1.5 + 228 * T ** -2.5) + \
                3.15e-29 * np.exp(-326 / T) * (0.5 * T ** -1.5 + 326 * T ** -2.5)
        dH2dT = 8.643e-15 * np.exp(-118400 / T) * T ** -2
        dCII_mdT = 2.755e-22 * np.exp(-6.2e4 / T) * (-0.5 * T ** -1.5 + 6.2e4 * T ** -2.5) + \
                   7.428e-23 * np.exp(-6.2e4 / T) * T ** -2
        dOI_mdT = 1.274e-27 * np.exp(-2.3e4 / T) * (-0.57 * T ** -0.43 + 2.3e4 * T ** -1.43) + \
                  1.095e-23 * np.exp(-2.3e4 / T) * T ** -2 + \
                  1.331e-27 * np.exp(-4.9e4 / T) * (-0.57 * T ** -0.43 + 4.9e4 * T ** -1.43) + \
                  4.969e-23 * np.exp(-4.9e4 / T) * T ** -2 + \
                  4.406e-27 * np.exp(-2.6e4 / T) * (-0.57 * T ** -0.43 + 2.6e4 * T ** -1.43) + \
                  4.199e-23 * np.exp(-2.6e4 / T) * T ** -2
        dRecdT = np.log(1.704 / T ** 0.068) * ((G0 * T ** 0.5) / nH) ** power / T ** 0.188 * \
                 (-2.340e-32 * T ** 0.06 * np.log(G0 * T ** 0.5 / nH) + 4.371e-31 * T ** 0.128 + 1.182e-31 * T ** 0.06)
    elif x <= 1.4e-4:
        dCIIdT = 1.394e-7 * np.exp(-92 / T) * (-0.5 * T ** -1.5 + 92 * T ** -2.5)
        dOIdT = 1.08e-28 * np.exp(-228 / T) * (0.5 * T ** -1.5 + 228 * T ** -2.5) + \
                3.15e-29 * np.exp(-326 / T) * (0.5 * T ** -1.5 + 326 * T ** -2.5)
        dH2dT = 1.21e-27 * np.exp(-118400 / T) * T ** -2
        dCII_mdT = 3.857e-25 * np.exp(-6.2e4 / T) * (-0.5 * T ** -1.5 + 6.2e4 * T ** -2.5) + \
                   7.428e-23 * np.exp(-6.2e4 / T) * T ** -2
        dOI_mdT = 1.789e-30 * np.exp(-2.3e4 / T) * (-0.57 * T ** -0.43 + 2.3e4 * T ** -1.43) + \
                  1.095e-23 * np.exp(-2.3e4 / T) * T ** -2 + \
                  1.863e-30 * np.exp(-4.9e4 / T) * (-0.57 * T ** -0.43 + 4.9e4 * T ** -1.43) + \
                  4.969e-23 * np.exp(-4.9e4 / T) * T ** -2 + \
                  6.165e-30 * np.exp(-2.6e4 / T) * (-0.57 * T ** -0.43 + 2.6e4 * T ** -1.43) + \
                  4.199e-23 * np.exp(-2.6e4 / T) * T ** -2
        dRecdT = np.log(6.567 / T ** 0.068) * ((G0 * T ** 0.5) / nH) ** power / T ** 0.188 * \
                 (-3.276e-35 * T ** 0.06 * np.log(G0 * T ** 0.5 / nH) + 6.119e-34 * T ** 0.128 - 4.982e-35 * T ** 0.06)
    else:
        dCIIdT = -1.511e-26 * np.exp(-92 / T) * (-0.25 * T ** -1.25 + 92 * T ** -2.25) / nH + \
                 2.06e-27 * np.exp(-92 / T) * (0.07 * T ** -0.93 + 92 * T ** -1.93)
        dOIdT = 1.08e-28 * np.exp(-228 / T) * (0.5 * T ** -1.5 + 228 * T ** -2.5) + \
                3.15e-29 * np.exp(-326 / T) * (0.5 * T ** -1.5 + 326 * T ** -2.5)
        dH2dT = 1.1081e-21 * np.exp(-118400 / T) * (0.25 * T ** -0.75 + 118400 * T ** -1.75) / nH
        dCII_mdT = 4.182e-24 * np.exp(-6.2e4 / T) * (-0.25 * T ** -1.25 + 6.2e4 * T ** -2.25) / nH + \
                   7.426e-23 * np.exp(-6.2e4 / T) * T ** -2
        dOI_mdT = 1.913e-29 * np.exp(-2.3e4 / T) * (0.82 * T ** -0.18 + 2.3e4 * T ** -1.18) / nH + \
                  1.096e-23 * np.exp(-2.3e4 / T) * T ** -2 + \
                  2.020e-29 * np.exp(-4.9e4 / T) * (0.82 * T ** -0.18 + 4.9e4 * T ** -1.18) / nH + \
                  4.969e-23 * np.exp(-4.9e4 / T) * T ** -2 + \
                  6.687e-29 * np.exp(-2.6e4 / T) * (0.82 * T ** -0.18 + 2.6e4 * T ** -1.18) / nH + \
                  4.199e-23 * np.exp(-2.6e4 / T) * T ** -2
        dRecdT = np.log(4.803 / T ** 0.068) * (G0 * T ** 0.25) ** power / (nH * T ** 1.068) * \
                 (-3.552e-34 * T ** 1.19 * np.log(G0 * T ** 0.25) + 8.400e-33 * T ** 1.258 - 9.995e-34 * T ** 1.19)
    return dCIIdT + dOIdT + dH2dT + dCII_mdT + dOI_mdT + dRecdT


def dheatingdrho(nH, T, x):
    G0 = 1 / 1.7
    if x >= 0.1:
        result = G0 * ((7.512e-28 * G0 ** 0.73 * T ** 0.365) / (
                nH ** 1.73 * ((0.021 * G0 ** 0.73 * T ** 0.365) / nH ** 0.73 + 1) ** 2)
                       + (1.173e-31 * G0 * T ** 1.2) / (nH ** 2 * ((0.002 * G0 * T ** 0.5) / nH + 1) ** 2))
    elif x <= 1.4e-4:
        result = G0 * ((9.315e-26 * G0 ** 0.73 * T ** 0.365) / (
                nH ** 1.73 * ((2.604 * G0 ** 0.73 * T ** 0.365) / nH ** 0.73 + 1) ** 2)
                       + (8.327e-29 * G0 * T ** 1.2) / (nH ** 2 * ((1.42 * G0 * T ** 0.5) / nH + 1) ** 2))
    else:
        result = -(1.635e-26 * G0 ** 1.73 * T ** 0.1825) / \
                 (nH ** 0.27 * (0.457 * G0 ** 0.73 * nH ** 0.73 * T ** 0.1825 + 1) ** 2)

    return result


def dheatingdT(nH, T, x):
    G0 = 1 / 1.7
    if x >= 0.1:
        result = G0 * (-(3.756e-28 * G0 ** 0.73) / (
                nH ** 0.73 * T ** 0.635 * ((0.021 * G0 ** 0.73 * T ** 0.365) / nH ** 0.73 + 1) ** 2)
                       + 4.105e-29 / (T ** 0.3 * ((0.002 * G0 * T ** 0.5) / nH + 1))
                       - (5.864e-32 * G0 * T ** 0.2) / (nH * ((0.002 * G0 * T ** 0.5) / nH + 1) ** 2))
    elif x <= 1.4e-4:
        result = G0 * (-(4.657e-26 * G0 ** 0.73) / (
                nH ** 0.73 * T ** 0.635 * ((2.604 * G0 ** 0.73 * T ** 0.365) / nH ** 0.73 + 1) ** 2)
                       + 4.105e-29 / (T ** 0.3 * ((1.42 * G0 * T ** 0.5) / nH + 1))
                       - (4.163e-29 * G0 * T ** 0.2) / (nH * ((1.42 * G0 * T ** 0.5) / nH + 1) ** 2))
    else:
        result = G0 * (-(4.087e-27 * G0 ** 0.73 * nH ** 0.73) / (
                T ** 0.8175 * (0.457 * G0 ** 0.73 * nH ** 0.73 * T ** 0.1825 + 1) ** 2)
                       + 4.105e-29 / (T ** 0.3 * (0.132 * G0 * T ** 0.5 + 1))
                       - (3.870e-30 * G0 * T ** 0.2) / (0.132 * G0 * T ** 0.5 + 1) ** 2)

    return result


def calculate_Gamma(m, n):
    gamma = 5 / 3
    c = 2 * (3 - 2 * m) / (5 - 2 * m - 2 * n)
    d = ((3 - 2 * m) / (5 - 2 * m - 2 * n) - 1) * 2
    Gamma1 = (gamma * 2 * (d / 3 + 1) + (n - m) * (c - (gamma - 1) * d)) / (
            2 * (d / 3 + 1) + (n - 1) * (c - (gamma - 1) * d))
    Gamma2 = 1 + (3 - 2 * m) / (2 * (n - 1))
    Gamma3 = 1 + (2 - m) / (n - 1)
    return Gamma1, Gamma2, Gamma3


def plot_m(fig, ax, nH_list, T_list, m_matrix):
    print("Extrema of m")
    print(np.min(m_matrix))
    print(np.max(m_matrix))
    print()
    m_matrix = np.flip(m_matrix.T, axis=0)
    im = ax.imshow(m_matrix, origin='lower', cmap="rainbow", vmin=0.1, vmax=10,
                   extent=(nH_list[0], nH_list[-1], T_list[0], T_list[-1]))
    ax.set_title("m(nH, T)")
    ax.set_xlabel(r'$\rm{log}nH \rm{cm^{-3}}$')
    ax.set_ylabel(r'$\rm{log}T \rm{K}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.colorbar(im, ax=ax, extend='both', label='m')
    return


def plot_n(fig, ax, nH_list, T_list, n_matrix):
    print("Extrema of n")
    print(np.min(n_matrix))
    print(np.max(n_matrix))
    print()
    n_matrix = np.flip(n_matrix.T, axis=0)
    im = ax.imshow(n_matrix, origin='lower', cmap="rainbow", vmin=0.01, vmax=0.1,
                   extent=(nH_list[0], nH_list[-1], T_list[0], T_list[-1]))
    ax.set_title("n(nH, T)")
    ax.set_xlabel(r'$\rm{log}nH\;\rm{cm^{-3}}$')
    ax.set_ylabel(r'$\rm{log}T\;\rm{K}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.colorbar(im, ax=ax, extend='both', label='n')
    return


def plot_Gamma1(fig, ax, nH_list, T_list, Gamma1):
    print("Extrema of Gamma_1")
    print(np.min(Gamma1))
    print(np.max(Gamma1))
    print()
    Gamma1 = np.flip(Gamma1.T, axis=0)
    im = ax.imshow(Gamma1, origin='lower', cmap="rainbow", vmin=0.5, vmax=5,
                   extent=(nH_list[0], nH_list[-1], T_list[0], T_list[-1]))
    ax.set_title(r"$\Gamma_1(nH,T)$")
    ax.set_xlabel(r"$\rm{log}nH\;\rm{cm^{-3}}$")
    ax.set_ylabel(r"$\rm{log}T\;\rm{K}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_1$")
    return


def plot_Gamma2(fig, ax, nH_list, T_list, Gamma2):
    print("Extrema of Gamma2")
    print(np.min(Gamma2))
    print(np.max(Gamma2))
    print()
    Gamma2 = np.flip(Gamma2.T, axis=0)
    im = ax.imshow(Gamma2, origin="lower", cmap="rainbow", vmin=0.5, vmax=5,
                   extent=(nH_list[0], nH_list[-1], T_list[0], T_list[-1]))
    ax.set_title(r"$\Gamma_2(nH,T)$")
    ax.set_xlabel(r"$\rm{log}nH\;\rm{cm^{-3}}$")
    ax.set_ylabel(r"$\rm{log}T\;\rm{K}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_2$")
    return


def plot_Gamma3(fig, ax, nH_list, T_list, Gamma3):
    print("Extrema of Gamma_3")
    print(np.min(Gamma3))
    print(np.max(Gamma3))
    print()
    Gamma3 = np.flip(Gamma3.T, axis=0)
    im = ax.imshow(Gamma3, origin="lower", cmap="rainbow", vmin=0.5, vmax=5,
                   extent=(nH_list[0], nH_list[-1], T_list[0], T_list[-1]))
    ax.set_title(r"$\Gamma_3(nH,T)$")
    ax.set_xlabel(r"$\rm{log}nH\;\rm{cm^{-3}}$")
    ax.set_ylabel(r"$\rm{log}T\;\rm{K}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_3$")
    return


def plot_Gamma1_m1(fig, ax, nH_list, T_list, Gamma1_m1):
    return


def plot_Gamma2_m1(fig, ax, nH_list, T_list, Gamma2_m1):
    return


def plot_Gamma3_m1(fig, ax, nH_list, T_list, Gamma3_m1):
    return


def make_plot(nH_list, T_list, m_matrix, n_matrix, Gamma1, Gamma2, Gamma3):
    fig, ((ax0, ax1, ax2), (ax3, ax4, ax5)) = plt.subplots(nrows=2, ncols=3, figsize=(21, 12))

    plot_m(fig, ax0, nH_list, T_list, m_matrix)
    plot_n(fig, ax1, nH_list, T_list, n_matrix)
    plot_Gamma1(fig, ax3, nH_list, T_list, Gamma1)
    plot_Gamma2(fig, ax4, nH_list, T_list, Gamma2)
    plot_Gamma3(fig, ax5, nH_list, T_list, Gamma3)

    ## quiver plot
    # plot_Gamma1_m1(fig, ax3, nH_list, T_list, Gamma1 - 1)
    # plot_Gamma2_m1(fig, ax4, nH_list, T_list, Gamma2 - 1)
    # plot_Gamma3_m1(fig, ax5, nH_list, T_list, Gamma3 - 1)

    plt.tight_layout(pad=2.5)
    plt.show()
    # plt.savfig("IC_Analyze.png", dpi=600)
    return


'''
Function	: main
Description	: calculate m, n from ramses cooling and heating cooling
Note		: 1. x axis of m_matrix is temperature which from low to high
              2. y axis of m_matrix is number density which from low to high 
			  3. n_matrix, Gamma1, Gamma2, Gamma3 are same as above
Parameter	: None
Return		: None
'''


def main():
    nH_start = -2
    nH_end = 4
    nH_count = 100
    T_start = -2
    T_end = 4
    T_count = 100

    nH_list = np.logspace(nH_start, nH_end, nH_count)
    T_list = np.logspace(T_start, T_end, T_count)
    m_matrix = np.zeros(shape=(nH_count, T_count))
    n_matrix = np.zeros(shape=(nH_count, T_count))

    for i in range(0, nH_count):
        for j in range(0, T_count):
            nH = nH_list[i]
            T = T_list[j]
            ne = 2.4e-3 * (T / 100) ** 0.25 / 0.5
            x = ne / nH

            total_cooling_rate = np.abs(nH * heating(nH, T, x) - nH ** 2 * cooling(nH, T, x))
            m_matrix[i][j] = nH / total_cooling_rate * np.abs(
                heating(nH, T, x) + nH * dheatingdrho(nH, T, x) - 2 * nH * cooling(nH, T, x) - nH ** 2 * dcoolingdrho(
                    nH, T, x))
            n_matrix[i][j] = T / total_cooling_rate * np.abs(nH * dheatingdT(nH, T, x) - nH ** 2 * dheatingdT(nH, T, x))

    Gamma1, Gamma2, Gamma3 = calculate_Gamma(m_matrix, n_matrix)
    make_plot(nH_list, T_list, m_matrix, n_matrix, Gamma1, Gamma2, Gamma3)


if __name__ == '__main__':
    main()
