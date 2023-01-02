import matplotlib.pyplot as plt

from module_calculate import *
from module_plot import *


def main():
    nH_start = -2
    nH_end = 10
    nH_count = 1000
    nH_array = np.logspace(nH_start, nH_end, nH_count)
    T_array = [10, 30, 50, 70, 100, 300, 500, 700, 1000, 3000, 5000, 7000]

    CII_cooling_rate = calculate_CII_cooling_rate(nH_array, T_array)
    OI_cooling_rate = calculate_OI_cooling_rate(nH_array, T_array)
    H2_cooling_rate = calculate_H2_cooling_rate(nH_array, T_array)
    CII_line_cooling_rate = calculate_CII_line_cooling_rate(nH_array, T_array)
    OI_line_cooling_rate = calculate_OI_line_cooling_rate(nH_array, T_array)
    Recombination_cooling_rate = calculate_Recombination_cooling_rate(nH_array, T_array)
    UV_heating_rate = calculate_UV_heating_rate(nH_array, T_array)

    plt.figure(figsize=(20, 15))
    for i in range(0, len(T_array)):
        ax = plt.subplot(3, 4, i + 1)
        T = T_array[i]
        plot_individual_cooling_rate_vs_nH(ax, nH_array, T, CII_cooling_rate.T[i], OI_cooling_rate.T[i],
                                           H2_cooling_rate.T[i], CII_line_cooling_rate.T[i],
                                           OI_line_cooling_rate.T[i],
                                           Recombination_cooling_rate.T[i], UV_heating_rate.T[i])

    plt.tight_layout()
    plt.savefig('cooling_rate_in_diff_temp.pdf', dpi=600)
    plt.show()

    return


if __name__ == "__main__":
    main()
