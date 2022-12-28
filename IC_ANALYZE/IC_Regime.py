import matplotlib.pyplot as plt

from module_calculate import *
from module_plot import *


def plot_slope_line(ax, slope, intercepts):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_vals = np.array(xlim)
    for intercept in intercepts:
        y_vals = intercept + slope * x_vals
        ax.plot(x_vals, y_vals, linestyle=':', color='black')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def main():
    nH_start = -2
    nH_end = 10
    nH_count = 1000
    T_start = 0
    T_end = 4
    T_count = 500

    nH_array = np.logspace(nH_start, nH_end, nH_count)
    T_array = np.logspace(T_start, T_end, T_count)

    balance_temperature = calculate_equilibrium_temperature(nH_array, T_array)

    m_matrix = calculate_m(nH_array, T_array, nH_count, T_count)
    n_matrix = calculate_n(nH_array, T_array, nH_count, T_count)
    Gamma1, Gamma2, Gamma3 = calculate_Gamma(m_matrix, n_matrix)

    fig, ((ax0, ax1, ax2), (ax3, ax4, ax5)) = plt.subplots(nrows=2, ncols=3, figsize=(21, 12))
    extent = [np.log10(nH_array[0]), np.log10(nH_array[-1]), np.log10(T_array[0]), np.log10(T_array[-1])]
    nH_mesh, T_mesh = np.meshgrid(np.log10(nH_array), np.log10(T_array))

    plot_m(fig, ax0, extent, m_matrix)
    plot_m_contour(fig, ax0, nH_mesh, T_mesh, m_matrix)
    plot_equilibrium_temperature(fig, ax0, np.log10(nH_array), np.log10(balance_temperature))

    plot_n(fig, ax1, extent, n_matrix)
    plot_n_contour(fig, ax1, nH_mesh, T_mesh, n_matrix)
    plot_equilibrium_temperature(fig, ax1, np.log10(nH_array), np.log10(balance_temperature))

    plot_Gamma1(fig, ax3, extent, Gamma1)
    plot_Gamma1_m1(fig, ax3, nH_mesh, T_mesh, Gamma1 - 1)
    plot_equilibrium_temperature(fig, ax3, np.log10(nH_array), np.log10(balance_temperature))

    plot_Gamma2(fig, ax4, extent, Gamma2)
    plot_Gamma2_m1(fig, ax4, nH_mesh, T_mesh, Gamma2 - 1)
    plot_Gamma2_contour(fig, ax4, nH_mesh, T_mesh, Gamma2)
    plot_slope_line(ax4, 2 / 3, np.arange(-6, 4, 2))
    # plot_slope_line(ax4, 1 / 3, np.arange(-6, 4, 1))
    # plot_slope_line(ax4, 0, np.arange(-6, 4, 1))
    plot_equilibrium_temperature(fig, ax4, np.log10(nH_array), np.log10(balance_temperature))

    plot_Gamma3(fig, ax5, extent, Gamma3)
    plot_Gamma3_m1(fig, ax5, nH_mesh, T_mesh, Gamma3 - 1)
    plot_equilibrium_temperature(fig, ax5, np.log10(nH_array), np.log10(balance_temperature))

    plt.tight_layout(pad=2.5)
    plt.savefig("IC_Analyze.pdf", dpi=600)
    plt.show()


if __name__ == "__main__":
    main()
