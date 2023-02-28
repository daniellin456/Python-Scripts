import matplotlib.pyplot as plt

from module_calculate import *
from module_plot import *


def main():
    T_start = 0
    T_end = 5
    T_count = 500
    T_array = np.logspace(T_start, T_end, T_count)
    nH_array = [5e-2, 1e-1, 5e-1, 1, 5, 10, 50, 1e2, 5e2, 1e3, 5e3, 1e4, 5e4,
                1e5, 5e5, 1e6]

    CII_cooling_rate = calculate_CII_cooling_rate(nH_array, T_array)
    OI_cooling_rate = calculate_OI_cooling_rate(nH_array, T_array)
    H2_cooling_rate = calculate_H2_cooling_rate(nH_array, T_array)
    CII_line_cooling_rate = calculate_CII_line_cooling_rate(nH_array, T_array)
    OI_line_cooling_rate = calculate_OI_line_cooling_rate(nH_array, T_array)
    Recombination_cooling_rate = calculate_Recombination_cooling_rate(nH_array, T_array)
    UV_heating_rate = calculate_UV_heating_rate(nH_array, T_array)

    fig, axs = plt.subplots(nrows=4, ncols=4, sharex=True, sharey=True, figsize=(16, 16))
    for i, ax in zip(range(0, len(nH_array)), axs.ravel()):
        nH = nH_array[i]
        plot_individual_cooling_rate_vs_T(ax, nH, T_array, CII_cooling_rate[i], OI_cooling_rate[i],
                                          H2_cooling_rate[i], CII_line_cooling_rate[i], OI_line_cooling_rate[i],
                                          Recombination_cooling_rate[i], UV_heating_rate[i])

        # ax.plot([0, 5], [-30, -25])
        # ax.plot([0, 5], [-30, -20])

    plt.suptitle(
        r"$\rm{log_{10}} \; n\Gamma, \rm{log_{10}} \; n^2 \Lambda$" + " vs " + "T" + " in different density",
        fontsize=18)
    fig.text(0.5, 0.01, r"$\rm{log_{10} \; T \; (K)}$", fontsize=14, ha='center')
    fig.text(0.01, 0.5, r"$\rm{log_{10}} \; n\Gamma, \rm{log_{10}} \; n^2 \Lambda \; (\rm{ergs \; cm^{-3}\; s^{-1}})$",
             va='center', rotation='vertical', fontsize=14)
    plt.tight_layout(rect=(0.02, 0.02, 0.99, 0.99))
    plt.savefig('cooling_rate_in_diff_density.pdf', dpi=600)
    plt.show()

    return


if __name__ == "__main__":
    main()
