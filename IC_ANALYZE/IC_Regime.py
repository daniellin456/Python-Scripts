import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import *
from Cooling_Function import *


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

    return balance_temperature


def calculate_m_n(nH_array, T_array, nH_count, T_count):
    m_matrix = np.zeros(shape=(nH_count, T_count))
    n_matrix = np.zeros(shape=(nH_count, T_count))
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
    return m_matrix, n_matrix


def calculate_Gamma(m, n):
    gamma = 5 / 3
    c = 2 * (3 - 2 * m) / (5 - 2 * m - 2 * n)
    d = ((3 - 2 * m) / (5 - 2 * m - 2 * n) - 1) * 2
    Gamma1 = (gamma * 2 * (d / 3 + 1) + (n - m) * (c - (gamma - 1) * d)) / (
            2 * (d / 3 + 1) + (n - 1) * (c - (gamma - 1) * d))
    Gamma2 = 1 + (3 - 2 * m) / (2 * (n - 1))
    Gamma3 = 1 + (2 - m) / (n - 1)
    return Gamma1, Gamma2, Gamma3


def plot_equilibrium_temperature(fig, ax, nH_array, balance_temperature):
    ax.plot(nH_array, balance_temperature, linewidth=3, linestyle="--", color="black")
    return


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
    im = ax.imshow(Gamma1.T, origin='lower', cmap="rainbow", norm=SymLogNorm(linthresh=1, vmin=-10, vmax=10), extent=extent)
    ax.set_title(r"$\Gamma_1\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_1$")
    return


def plot_Gamma2(fig, ax, extent, Gamma2):
    print("Extrema of Gamma2: min:%24.14e, max: %24.14e" % (np.min(Gamma2), np.max(Gamma2)))
    im = ax.imshow(Gamma2.T, origin="lower", cmap="rainbow", norm=SymLogNorm(linthresh=1, vmin=-10, vmax=10), extent=extent)
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


def plot_Gamma1_m1(fig, ax, nH_array, T_array, Gamma1_m1):
    nH_mesh, T_mesh = np.meshgrid(np.log10(nH_array), np.log10(T_array))
    ax.streamplot(nH_mesh, T_mesh, np.ones(shape=nH_mesh.shape), Gamma1_m1.T, color="white")
    return


def plot_Gamma2_m1(fig, ax, nH_array, T_array, Gamma2_m1):
    nH_mesh, T_mesh = np.meshgrid(np.log10(nH_array), np.log10(T_array))
    ax.streamplot(nH_mesh, T_mesh, np.ones(shape=nH_mesh.shape), Gamma2_m1.T, color="white")
    return


def plot_Gamma3_m1(fig, ax, nH_array, T_array, Gamma3_m1):
    nH_mesh, T_mesh = np.meshgrid(np.log10(nH_array), np.log10(T_array))
    ax.streamplot(nH_mesh, T_mesh, np.ones(shape=nH_mesh.shape), Gamma3_m1.T, color="white")
    return


def make_plot(nH_array, T_array, balance_temperature, m_matrix, n_matrix, Gamma1, Gamma2, Gamma3):
    fig, ((ax0, ax1, ax2), (ax3, ax4, ax5)) = plt.subplots(nrows=2, ncols=3, figsize=(21, 12))
    extent = [np.log10(nH_array[0]), np.log10(nH_array[-1]), np.log10(T_array[0]), np.log10(T_array[-1])]
    nH_mesh, T_mesh = np.meshgrid(nH_array, T_array)

    plot_m(fig, ax0, extent, m_matrix)
    plot_equilibrium_temperature(fig, ax0, np.log10(nH_array), np.log10(balance_temperature))

    plot_n(fig, ax1, extent, n_matrix)
    plot_equilibrium_temperature(fig, ax1, np.log10(nH_array), np.log10(balance_temperature))

    plot_Gamma1(fig, ax3, extent, Gamma1)
    plot_equilibrium_temperature(fig, ax3, np.log10(nH_array), np.log10(balance_temperature))

    plot_Gamma2(fig, ax4, extent, Gamma2)
    plot_equilibrium_temperature(fig, ax4, np.log10(nH_array), np.log10(balance_temperature))

    plot_Gamma3(fig, ax5, extent, Gamma3)
    plot_equilibrium_temperature(fig, ax5, np.log10(nH_array), np.log10(balance_temperature))

    ## quiver plot
    plot_Gamma1_m1(fig, ax3, nH_array, T_array, Gamma1 - 1)
    plot_Gamma2_m1(fig, ax4, nH_array, T_array, Gamma2 - 1)
    plot_Gamma3_m1(fig, ax5, nH_array, T_array, Gamma3 - 1)

    plt.tight_layout(pad=2.5)
    plt.show()
    # plt.savfig("IC_Analyze.png", dpi=600)
    return


def main():
    nH_start = -2
    nH_end = 6
    nH_count = 500
    T_start = 0
    T_end = 4
    T_count = 500

    nH_array = np.logspace(nH_start, nH_end, nH_count)
    T_array = np.logspace(T_start, T_end, T_count)
    balance_temperature = calculate_equilibrium_temperature(nH_array, T_array)
    m_matrix, n_matrix = calculate_m_n(nH_array, T_array, nH_count, T_count)
    Gamma1, Gamma2, Gamma3 = calculate_Gamma(m_matrix, n_matrix)
    make_plot(nH_array, T_array, balance_temperature, m_matrix, n_matrix, Gamma1, Gamma2, Gamma3)


if __name__ == "__main__":
    main()
