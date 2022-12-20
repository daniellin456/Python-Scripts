from module_plot import *
from module_calculate import *


def main():
    nH_start = -2
    nH_end = 6
    nH_count = 500
    nH_array = np.logspace(nH_start, nH_end, nH_count)
    T_array = np.asarray([10, 50, 100, 500, 1000])
    T_count = len(T_array)

    m_matrix = calculate_m(nH_array, T_array, nH_count, T_count)
    dLambda_dRho_by_np_diff = calculate_dLambda_dRho_by_np_diff(nH_array, T_array)
    dLambda_dRho_by_analytical = calculate_dLambda_dRho_by_analytical(nH_array, T_array)

    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, figsize=(24, 8))
    extent = [np.log10(nH_array[0]), np.log10(nH_array[-1]), np.log10(T_array[0]), np.log10(T_array[-1])]

    plot_m(fig, ax0, extent, m_matrix)
    plot_dLambdadRho_by_np_diff(fig, ax1, nH_array, T_array, dLambda_dRho_by_np_diff)
    plot_dLambdadRho_by_analytical(fig, ax2, nH_array, T_array, dLambda_dRho_by_analytical)
    plt.tight_layout()
    plt.show()

    return


if __name__ == "__main__":
    main()
