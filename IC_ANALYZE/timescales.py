import numpy as np
import matplotlib.pyplot as plt
from IC_ANALYZE.IC_Regime import *

kb = 1.3606e-16
gamma = 5/3
G = 6.67e-8

def calculate_cooling_timescale(nH, T, x):
    total_cooling_rate = cooling(nH, T, x)
    cooling_timescale = nH * kb * T / ( nH ** 2 * total_cooling_rate * (gamma - 1))
    return cooling_timescale


def calculate_heating_timescale(nH, T, x):
    total_heating_rate = heating(nH, T, x)
    heating_timescale = nH * kb * T / (nH * total_heating_rate * (gamma - 1))
    return heating_timescale

def calculate_free_fall_timescale(nH, T, x):
    free_fall_timescale = 1/np.sqrt(nH * G)

def main():
    nH_start = -2
    nH_end = 6
    nH_count = 100
    T_start = -2
    T_end = 6
    T_count = 100

    nH_array = np.logspace(nH_start, nH_end, nH_count)
    T_array = np.logspace(T_start, T_end, T_count)
    time = np.zeros(shape=(nH_count, T_count))

    day = 86400
    year = 365 * day
    kyr = 1e3 * year
    myr = 1e6 * year

    for i in range(0, nH_count):
        for j in range(0, T_count):
            nH = nH_array[i]
            T = T_array[j]
            ne = 2.4e-3 * ((T / 100.0) ** 0.25) / 0.5
            x = ne / nH

            heating_timescale = calculate_heating_timescale(nH, T, x)
            cooling_timescale = calculate_cooling_timescale(nH, T, x)
            time[i, j] = heating_timescale - cooling_timescale

    sc = plt.imshow(time.T/myr, origin='lower', extent=(nH_start, nH_end, T_start, T_end), norm=SymLogNorm(linthresh=0.01, vmin=-1e3, vmax=1e3), cmap="RdBu")
    plt.colorbar(sc)
    plt.show()



    return


if __name__ == "__main__":
    main()
