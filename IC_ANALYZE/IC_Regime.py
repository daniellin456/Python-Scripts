import matplotlib.pyplot as plt

from module_calculate import *
from module_plot import *


def intergrated_plot(nH_array, T_array, m_matrix, n_matrix, Gamma1, Gamma2, Gamma3, a, b, c, d, balance_temperature, total_cooling_power):
    gamma = 5 / 3
    fig, ((ax0, ax1, ax2), (ax3, ax4, ax5), (ax6, ax7, ax8)) = plt.subplots(nrows=3, ncols=3, figsize=(21, 18))
    extent = [np.log10(nH_array[0]), np.log10(nH_array[-1]), np.log10(T_array[0]), np.log10(T_array[-1])]
    nH_mesh, T_mesh = np.meshgrid(np.log10(nH_array), np.log10(T_array))

    plot_m(fig, ax0, extent, m_matrix)
    plot_m_contour(fig, ax0, nH_mesh, T_mesh, m_matrix)
    plot_equilibrium_temperature(fig, ax0, np.log10(nH_array), np.log10(balance_temperature))

    plot_n(fig, ax1, extent, n_matrix)
    plot_n_contour(fig, ax1, nH_mesh, T_mesh, n_matrix)
    plot_equilibrium_temperature(fig, ax1, np.log10(nH_array), np.log10(balance_temperature))

    plot_total_cooling_rate(fig, ax2, extent, total_cooling_power)
    plot_equilibrium_temperature(fig, ax2, np.log10(nH_array), np.log10(balance_temperature))

    plot_ax3(fig, ax3, extent, 3 - 2 * a)
    plot_ax3_contour(fig, ax3, nH_mesh, T_mesh, 3 - 2 * a)
    plot_ax4(fig, ax4, extent, c / (gamma - 1) - d)
    plot_ax4_contour(fig, ax4, nH_mesh, T_mesh, c / (gamma - 1) - d)
    plot_ax5(fig, ax5, extent, c / (gamma - 1) + b + 2)
    plot_ax5_contour(fig, ax5, nH_mesh, T_mesh, c / (gamma - 1) + b + 2)

    plot_Gamma1(fig, ax6, extent, Gamma1)
    plot_Gamma1_m1(fig, ax6, nH_mesh, T_mesh, Gamma1 - 1)
    plot_equilibrium_temperature(fig, ax6, np.log10(nH_array), np.log10(balance_temperature))

    plot_Gamma2(fig, ax7, extent, Gamma2)
    plot_Gamma2_m1(fig, ax7, nH_mesh, T_mesh, Gamma2 - 1)
    plot_Gamma2_contour(fig, ax7, nH_mesh, T_mesh, Gamma2)
    plot_slope_line(ax7, 2 / 3, np.arange(-6, 4, 2), 'black')
    plot_slope_line(ax7, 1 / 3, np.arange(-6, 4, 2), 'red')
    plot_equilibrium_temperature(fig, ax7, np.log10(nH_array), np.log10(balance_temperature))

    plot_Gamma3(fig, ax8, extent, Gamma3)
    plot_Gamma3_m1(fig, ax8, nH_mesh, T_mesh, Gamma3 - 1)
    plot_equilibrium_temperature(fig, ax8, np.log10(nH_array), np.log10(balance_temperature))

    plt.tight_layout(pad=2.5)
    plt.savefig("IC_Analyze.pdf", dpi=600)
    plt.show()
    return


def individual_plot(nH_array, T_array, m_matrix, n_matrix, Gamma1, Gamma2, Gamma3, a, b, c, d, balance_temperature, total_cooling_power):
    fig_x_length = 5
    fig_y_length = 3
    gamma = 5 / 3
    extent = [np.log10(nH_array[0]), np.log10(nH_array[-1]), np.log10(T_array[0]), np.log10(T_array[-1])]
    nH_mesh, T_mesh = np.meshgrid(np.log10(nH_array), np.log10(T_array))

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_x_length, fig_y_length))
    plot_total_cooling_rate(fig, ax, extent, total_cooling_power)
    plot_equilibrium_temperature(fig, ax, np.log10(nH_array), np.log10(balance_temperature))
    plt.tight_layout()
    fig.savefig('equilibrium_temperature.pdf', dpi=600)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_x_length, fig_y_length))
    plot_m(fig, ax, extent, m_matrix)
    plot_m_contour(fig, ax, nH_mesh, T_mesh, m_matrix)
    plot_equilibrium_temperature(fig, ax, np.log10(nH_array), np.log10(balance_temperature))
    plt.tight_layout()
    fig.savefig('m.pdf', dpi=600)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_x_length, fig_y_length))
    plot_n(fig, ax, extent, n_matrix)
    plot_n_contour(fig, ax, nH_mesh, T_mesh, n_matrix)
    plot_equilibrium_temperature(fig, ax, np.log10(nH_array), np.log10(balance_temperature))
    plt.tight_layout()
    fig.savefig('n.pdf', dpi=600)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_x_length, fig_y_length))
    plot_ax3(fig, ax, extent, 3 - 2 * a)
    plot_ax3_contour(fig, ax, nH_mesh, T_mesh, 3 - 2 * a)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_x_length, fig_y_length))
    plot_ax4(fig, ax, extent, c / (gamma - 1) - d)
    plot_ax4_contour(fig, ax, nH_mesh, T_mesh, c / (gamma - 1) - d)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_x_length, fig_y_length))
    plot_ax5(fig, ax, extent, c / (gamma - 1) + b + 2)
    plot_ax5_contour(fig, ax, nH_mesh, T_mesh, c / (gamma - 1) + b + 2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_x_length, fig_y_length))
    plot_Gamma1(fig, ax, extent, Gamma1)
    plot_Gamma1_m1(fig, ax, nH_mesh, T_mesh, Gamma1 - 1)
    plot_equilibrium_temperature(fig, ax, np.log10(nH_array), np.log10(balance_temperature))
    fig.savefig('Gamma1.pdf', dpi=600)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_x_length, fig_y_length))
    plot_Gamma2(fig, ax, extent, Gamma2)
    plot_Gamma2_m1(fig, ax, nH_mesh, T_mesh, Gamma2 - 1)
    plot_Gamma2_contour(fig, ax, nH_mesh, T_mesh, Gamma2)
    plot_slope_line(ax, 2 / 3, np.arange(-6, 4, 2), 'black')
    plot_slope_line(ax, 1 / 3, np.arange(-6, 4, 2), 'red')
    plot_equilibrium_temperature(fig, ax, np.log10(nH_array), np.log10(balance_temperature))
    fig.savefig('Gamma2.pdf', dpi=600)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_x_length, fig_y_length))
    plot_Gamma3(fig, ax, extent, Gamma3)
    plot_Gamma3_m1(fig, ax, nH_mesh, T_mesh, Gamma3 - 1)
    plot_equilibrium_temperature(fig, ax, np.log10(nH_array), np.log10(balance_temperature))
    fig.savefig('Gamma3.pdf', dpi=600)

    return


def main():
    nH_start = -2
    nH_end = 10
    nH_count = 500
    T_start = 0
    T_end = 4
    T_count = 500

    nH_array = np.logspace(nH_start, nH_end, nH_count)
    T_array = np.logspace(T_start, T_end, T_count)

    balance_temperature = calculate_equilibrium_temperature(nH_array, T_array)
    total_cooling_power = calculate_total_cooling_power(nH_array, T_array, nH_count, T_count)

    m_matrix = calculate_m(nH_array, T_array, nH_count, T_count)
    n_matrix = calculate_n(nH_array, T_array, nH_count, T_count)
    a = calculate_a(m_matrix, n_matrix)
    b = calculate_b(m_matrix, n_matrix)
    c = calculate_c(m_matrix, n_matrix)
    d = calculate_d(m_matrix, n_matrix)
    Gamma1 = calculate_Gamma1(m_matrix, n_matrix)
    Gamma2 = calculate_Gamma2(m_matrix, n_matrix)
    Gamma3 = calculate_Gamma3(m_matrix, n_matrix)

    intergrated_plot(nH_array, T_array, m_matrix, n_matrix, Gamma1, Gamma2, Gamma3, a, b, c, d, balance_temperature, total_cooling_power)
    individual_plot(nH_array, T_array, m_matrix, n_matrix, Gamma1, Gamma2, Gamma3, a, b, c, d, balance_temperature, total_cooling_power)
    return


if __name__ == "__main__":
    main()
