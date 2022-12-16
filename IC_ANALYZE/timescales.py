import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import *
from Cooling_Function import *
from IC_Regime import *

kb = 1.3606e-16
gamma = 5 / 3
G = 6.67e-8
mu = 2.31
mp = 1.66e-24


def main():
    nH_start = -2
    nH_end = 6
    nH_count = 100
    T_start = 0
    T_end = 4
    T_count = 100

    nH_array = np.logspace(nH_start, nH_end, nH_count)
    T_array = np.logspace(T_start, T_end, T_count)
    timescale_ratio = np.zeros(shape=(nH_count, T_count))

    for i in range(0, nH_count):
        for j in range(0, T_count):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH

            total_cooling_rate = np.abs(nH * heating(nH, T, x) - nH ** 2 * cooling(nH, T, x))
            cooling_heating_timescale = nH * kb * T / ((gamma - 1) * total_cooling_rate)
            free_fall_Timescale = np.sqrt(1 / (G * nH * mu * mp))
            timescale_ratio[i, j] = cooling_heating_timescale / free_fall_Timescale
    print(
        "Extrema of  timescale_ratio, min: %24.14e, max: %24.14e" % (np.min(timescale_ratio), np.max(timescale_ratio)))
    sc = plt.imshow(timescale_ratio.T, origin='lower', extent=(nH_start, nH_end, T_start, T_end), cmap="RdBu",
                    norm=LogNorm(vmin=1e-3, vmax=1e3))
    plt.plot(np.log10(nH_array), np.log10(calculate_equilibrium_temperature(nH_array, T_array)), color="black", linestyle="--")
    plt.colorbar(sc)
    plt.show()

    return


if __name__ == "__main__":
    main()
