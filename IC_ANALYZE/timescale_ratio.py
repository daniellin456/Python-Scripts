import matplotlib.pyplot as plt

from IC_Regime import *
from module_style import *

kb = 1.3606e-16
gamma = 5 / 3
G = 6.67e-8
mu = 2.31
mp = 1.66e-24


def main():
    plot_style()

    nH_start = -2
    nH_end = 10
    nH_count = 500
    T_start = 0
    T_end = 4
    T_count = 500

    nH_array = np.logspace(nH_start, nH_end, nH_count)
    T_array = np.logspace(T_start, T_end, T_count)
    timescale_ratio = np.zeros(shape=(nH_count, T_count))

    nH_mesh, T_mesh = np.meshgrid(np.log10(nH_array), np.log10(T_array))
    for i in range(0, nH_count):
        for j in range(0, T_count):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH

            total_cooling_rate = np.abs(nH * heating(nH, T, x) - nH ** 2 * cooling(nH, T, x))
            cooling_heating_timescale = nH * kb * T / ((gamma - 1) * total_cooling_rate)
            free_fall_Timescale = np.sqrt(3 * np.pi / (32 * G * nH * mu * mp))
            timescale_ratio[i, j] = cooling_heating_timescale / free_fall_Timescale
    print(
        "Extrema of  timescale_ratio, min: %24.14e, max: %24.14e" % (np.min(timescale_ratio), np.max(timescale_ratio)))
    # fig, ax = plt.figure(figsize=(6, 4))
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))
    sc = plt.imshow(timescale_ratio.T, origin='lower', extent=(nH_start, nH_end, T_start, T_end), cmap="RdBu",
                    norm=LogNorm(vmin=1e-3, vmax=1e3))
    plt.plot(np.log10(nH_array), np.log10(calculate_equilibrium_temperature(nH_array, T_array)), color="black",
             linestyle="--")
    CS = plt.contour(nH_mesh, T_mesh, timescale_ratio.T, levels=[0.001, 0.01, 0.1], colors=('cyan', 'yellow', 'green'))
    plt.clabel(CS, CS.levels, inline=True, fontsize=14)
    plt.title(r"$t_{\rm{cool}} / t_{\rm{ff}}$")
    plt.xlabel(r"$\rm{log_{10}} n_H \; (\rm{cm^{-3}})$")
    plt.ylabel(r"$\rm{log_{10}} T \; (\rm{K})$")
    plt.xticks([-2, 0, 2, 4, 6, 8, 10])
    plt.yticks([0, 1, 2, 3, 4])
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    cb = plt.colorbar(sc, location="bottom", pad=0.3, extend="both", ticks=[1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3])
    cb.ax.tick_params(labelsize=14)
    cb.ax.minorticks_on()
    plt.savefig("timescale_ratio.pdf")
    plt.tight_layout()
    plt.show()

    return


if __name__ == "__main__":
    main()
