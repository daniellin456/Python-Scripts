from module_calculate import *
from module_plot import *


def main():
    nH_start = -2
    nH_end = 6
    nH_count = 1000
    nH_array = np.logspace(nH_start, nH_end, nH_count)
    T_array = np.asarray([10, 50, 100, 500, 1000])

    dLambda_dRho_by_np_diff = calculate_dLambda_dRho_by_np_diff(nH_array, T_array)
    dLambda_dRho_by_analytical = calculate_dLambda_dRho_by_analytical(nH_array, T_array)
    dGamma_dRho_by_np_diff = calculate_dGamma_dRho_by_np_diff(nH_array, T_array)
    dGamma_dRho_by_np_analytical = calculate_dGamma_dRho_by_analytical(nH_array, T_array)
    dTotal_Lambda_dRho_by_np_diff = calculate_dTotal_Lambda_dRho_by_np_diff(nH_array, T_array)
    dTotal_Lambda_dRho_by_analytical = calculate_dTotal_Lambda_dRho_by_analytical(nH_array, T_array)
    ionization = calculate_ionization(nH_array, T_array)

    fig, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(nrows=3, ncols=2, figsize=(10, 15))

    plot_dTotal_Lambda_dRho_by_np_diff(fig, ax0, nH_array, T_array, dTotal_Lambda_dRho_by_np_diff)
    plot_dTotal_Lambda_dRho_by_analytical(fig, ax1, nH_array, T_array, dTotal_Lambda_dRho_by_analytical)

    plot_dLambda_dRho_by_np_diff(fig, ax2, nH_array, T_array, dLambda_dRho_by_np_diff)
    plot_ionization(fig, ax2, nH_array, T_array, ionization)

    plot_dLambda_dRho_by_analytical(fig, ax3, nH_array, T_array, dLambda_dRho_by_analytical)
    plot_ionization(fig, ax3, nH_array, T_array, ionization)

    plot_dGamma_dRho_by_np_diff(fig, ax4, nH_array, T_array, dGamma_dRho_by_np_diff)
    plot_ionization(fig, ax4, nH_array, T_array, ionization)

    plot_dGamma_dRho_by_analytical(fig, ax5, nH_array, T_array, dGamma_dRho_by_np_analytical)
    plot_ionization(fig, ax5, nH_array, T_array, ionization)

    plt.tight_layout()
    plt.show()

    return


if __name__ == "__main__":
    main()
