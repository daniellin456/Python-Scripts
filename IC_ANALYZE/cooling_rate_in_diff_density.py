import matplotlib.pyplot as plt

from module_calculate import *
from module_plot import *


def main():
    T_start = 0
    T_end = 5
    T_count = 500
    T_array = np.logspace(T_start, T_end, T_count)
    nH_array = [1e-2, 5e-2, 1e-1, 5e-1, 1, 5, 10, 50, 1e2, 5e2, 1e3, 5e3, 1e4, 5e4,
                1e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9, 5e9]

    CII_cooling_rate = calculate_CII_cooling_rate(nH_array, T_array)
    OI_cooling_rate = calculate_OI_cooling_rate(nH_array, T_array)
    H2_cooling_rate = calculate_H2_cooling_rate(nH_array, T_array)
    CII_line_cooling_rate = calculate_CII_line_cooling_rate(nH_array, T_array)
    OI_line_cooling_rate = calculate_OI_line_cooling_rate(nH_array, T_array)
    Recombination_cooling_rate = calculate_Recombination_cooling_rate(nH_array, T_array)
    UV_heating_rate = calculate_UV_heating_rate(nH_array, T_array)

    plt.figure(figsize=(20, 30))
    for i in range(0, len(nH_array)):
        ax = plt.subplot(6, 4, i + 1)
        nH = nH_array[i]
        plot_individual_cooling_rate_vs_T(ax, nH, T_array, CII_cooling_rate[i], OI_cooling_rate[i],
                                          H2_cooling_rate[i], CII_line_cooling_rate[i], OI_line_cooling_rate[i],
                                          Recombination_cooling_rate[i], UV_heating_rate[i])

        # ax.plot([0, 5], [-30, -25])
        # ax.plot([0, 5], [-30, -20])

    plt.tight_layout()
    plt.savefig('cooling_rate_in_diff_density.pdf', dpi=600)
    plt.show()

    return


if __name__ == "__main__":
    main()
