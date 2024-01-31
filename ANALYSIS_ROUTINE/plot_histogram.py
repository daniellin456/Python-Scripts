import os
import argparse
import matplotlib
matplotlib.use('Agg')

import pymses
from pymses.sources.ramses.filename_utils import  search_valid_outputs
from pymses.filters import CellsToPoints
from scipy.optimize import fsolve
from module_histogram import *
from module_units import *
from module_plot import *
from module_style import *

G = 6.67e-8
mu = 2.31
mp = 1.67e-24
gamma = 1.6666667
kb = 1.3806e-16

def find_free_fall_time(target_folder_path):

    peak_density = 0.0
    if "SC01" in target_folder_path:
        peak_density = 1280.01295051676

    elif "SC0068" in target_folder_path:
        peak_density = 12946.7537409113

    elif "SC0047" in target_folder_path:
        peak_density = 118748.249427123

    free_fall_time = np.sqrt(3.0 * np.pi / 32.0 / 6.67e-8 / 2.31 / 1.66e-24 / peak_density)

    return free_fall_time

def find_eta(target_folder_path):
    eta = ""
    if "0.1_EMPIRICAL_COOLING" in target_folder_path:
        index = target_folder_path.index("0.1_EMPIRICAL_COOLING")
        eta = target_folder_path[index:index + len("0.1")]
    elif "0.01_EMPIRICAL_COOLING" in target_folder_path:
        index = target_folder_path.index("0.01_EMPIRICAL_COOLING")
        eta = target_folder_path[index:index + len("0.01")]
    elif "0.001_EMPIRICAL_COOLING" in target_folder_path:
        index = target_folder_path.index("0.001_EMPIRICAL_COOLING")
        eta = target_folder_path[index:index + len("0.001")]
    else:
        eta = "1.0"
    return eta

def find_ffsct(target_folder_path):

    ffsct = ""
    if "SC0.100" in target_folder_path:
        index = target_folder_path.index("0.100")
        ffsct = target_folder_path[index:index + len("0.100")]
    elif "0.068" in target_folder_path:
        index = target_folder_path.index("0.068")
        ffsct = target_folder_path[index:index + len("0.068")]
    elif "0.047" in target_folder_path:
        index = target_folder_path.index("0.047")
        ffsct = target_folder_path[index:index + len("0.047")]

    return ffsct

def find_T0(target_folder_path):
    T0 = ""
    if "10K" in target_folder_path:
        index = target_folder_path.index("10K")
        T0 = target_folder_path[index:index + len("10K")]
    elif "80K" in target_folder_path:
        index = target_folder_path.index("80K")
        T0 = target_folder_path[index:index + len("80K")]
    elif "400K" in target_folder_path:
        index = target_folder_path.index("400K")
        T0 = target_folder_path[index:index + len("400K")]

    return T0

def get_parametric_cooling_m_n(target_folder_path):

    m = ""
    n = ""
    lambda_0 = ""
    if "M" and "N" in target_folder_path:
        with open("/data/daniellin/RAMSES_RUNS/" + target_folder_path +"/run.nml", "r") as file:
            for line in file:
                if "power_m" in line:
                    index = line.index("power_m")
                    m = line[index+len("power_m")+1:index+len("power_m")+4]

                if "power_n" in line:
                    index = line.index("power_n")
                    n = line[index+len("power_n")+1:index+len("power_n")+4]

                if "lambda_0" in line:
                    index = line.index("lambda_0")
                    lambda_0 = line[index+len("lambda_0")+1:index+len("lambda_0")+10]
    return m, n, lambda_0

def get_slope(x1, y1, x2, y2):
    return (y2 - y1) / (x2 - x1)

def func(T, *args):
    nH, m, n, lambda_0 = args
    return np.sqrt(3*np.pi / 32.0 / G/ mu / mp / nH) - nH * kb * T / (gamma - 1) / lambda_0 / nH ** m / T ** n

def plot_t_ff_eq_t_cooling(ax, m, n, lambda_0):
    print(m, n, lambda_0)
    nH_start = -2
    nH_end = 10
    nH_count = 1000
    nH_array = np.logspace(nH_start, nH_end, nH_count)
    root = []

    for nH in nH_array:
        result = fsolve(func, 1.0, args=(nH, m, n, lambda_0))[0]
        root.append(result)

    slope = get_slope(np.log10(nH_array[0]), np.log10(root[0]), np.log10(nH_array[-1]), np.log10(root[-1]))
    ax.plot(np.log10(nH_array), np.log10(root), '--', lw=2, color='black', label=r"$\gamma = $" + str(round(slope+1, 3)))
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    return slope


def main(target_folder_paths, writedir, timestamps):

    for target_folder_path in target_folder_paths: # Loop FFSCT = 0.1, 0.068, 0.047

        free_fall_time = find_free_fall_time(target_folder_path)
        # ffsct = find_ffsct(target_folder_path)
        # T0 = find_T0(target_folder_path)
        # eta = find_eta(target_folder_path)
        m, n, lambda_0 = get_parametric_cooling_m_n(target_folder_path.split('/')[-1])

        # for timestamp in timestamps[target_folder_paths.index(target_folder_path)]: # Draw 3 different timestamp in each FFSCT
        for timestamp in timestamps:  # Draw 3 different timestamp in each FFSCT
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4.5))

            data = pymses.RamsesOutput(target_folder_path, timestamp)
            units = set_units(data)
            time_sec = data.info['time'] * units['unit_T_s']
            time_kyr = data.info['time'] * units['unit_T_Kyr']
            amr_data = data.amr_source(["rho", "T", "vel"])
            cell = CellsToPoints(amr_data).flatten()
            volume = cell.get_sizes() ** 3
            density = cell["rho"]
            temperature = cell["T"]
            mass = density * volume

            density_temperature_2d_histogram, density_edges, temperature_edges = calculate_density_temperature_2d_histogram(
                density, temperature, mass)
            nH, equilibrium_temperature = calculate_equilibrium_temperature()
            plot_density_temperature_2d_histogram_imshow(fig, ax, density_temperature_2d_histogram.T *
                                                         units["unit_M_Msun"],
                                                         density_edges + np.log10(units["unit_D_Hcc"]),
                                                         temperature_edges)

            if m != "" and n != "" and lambda_0 != "":
                # plot_t_ff_eq_t_cooling(ax, float(m), float(n), float(lambda_0))
                if m == "1.3" and n == "2.3":
                    plot_t_ff_eq_t_cooling(ax, 1.3, 2.3, 3.987e-32)
                elif m == "1.5" and n == "3.0":
                    plot_t_ff_eq_t_cooling(ax, 1.5, 3.0, 1.902e-33)
                elif m == "2.0" and n == "4.0":
                    plot_t_ff_eq_t_cooling(ax, 2.0, 4.0, 5.316e-36)

                # ax.set_title(r"$t_{ff}/t_{sc}$ = " + ffsct + ", m = " + m + ", n = " + n+ r", $T_0$ = " + T0)#, fontsize=8)
                ax.legend(loc="lower right", fontsize=10)

            else:
                plot_theoretical_equilibrium_temperature(ax, np.log10(nH), np.log10(equilibrium_temperature))
                # ax.set_title(r"$t_{ff}/t_{sc}$ = " + ffsct + r", $\eta$ = " + eta + r", $T_0$ = " + T0)


            # fig.suptitle("t = " + str(round(time_sec/free_fall_time, 3)) + r"$ t_{ff}$, " + str(round(time_kyr, 3)) + " kyr", fontsize=14, y=0.90)
            ax.text(9.5, 3.5, "t = " + str(round(time_sec/free_fall_time, 3)) + r"$ t_{ff}$, " + str(round(time_kyr, 3)) + " kyr", ha='right', va='center', fontsize=18)
            fig.tight_layout()
            time = str(round(time_kyr, 3))

            # if m != "" and n != "":
                # time = str(round(time_kyr, 3))
                # filename = "m_" + m + "_n_" + n + "_T0_" + T0 + "_ffsct_" + ffsct + "_t_" + time + "_Hist.pdf"
                # fig.savefig(filename)
                # fig.savefig(writedir + "/Hist_t_" + time +".pdf")
            # else:
                # time = str(round(time_kyr, 3))
                # filename = "eta_" + eta + "_T0_" + T0 + "_ffsct_" + ffsct + "_t_" + time +"_Hist.pdf"
                # fig.savefig(filename)
                # fig.savefig(writedir + "/Hist_t_" + time + ".pdf")
            fig.savefig(writedir + "/Hist_" + format(timestamp, '05') + ".pdf")

    return

if __name__ == "__main__":

    plot_style()
    parser = argparse.ArgumentParser(description="plot profile")
    parser.add_argument("read_dir", type=str, help="RMASES output repository")
    parser.add_argument("-t", "--time", help="timestamp")
    args = parser.parse_args()

    readdir = [args.read_dir]
    writedir = "./figure/" + args.read_dir.split('/')[-1]
    timestamps = eval(args.time)

    # ROOT_PATH = "/data/daniellin/RAMSES_RUNS/"
    # TARGET_FOLDER_PATHS = [
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.100",
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.068",
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.047",
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_80K/FFSCT_0.100",
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_80K/FFSCT_0.068",
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_80K/FFSCT_0.047",
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_400K/FFSCT_0.100",
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_400K/FFSCT_0.068",
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_400K/FFSCT_0.047",
        # ROOT_PATH + "0.1_EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.100",
        # ROOT_PATH + "0.1_EMPIRICAL_COOLING/BG_TEMP_80K/FFSCT_0.100",
        # ROOT_PATH + "0.1_EMPIRICAL_COOLING/BG_TEMP_400K/FFSCT_0.100",
        # ROOT_PATH + "0.01_EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.100",
        # ROOT_PATH + "0.01_EMPIRICAL_COOLING/BG_TEMP_80K/FFSCT_0.100",
        # ROOT_PATH + "0.01_EMPIRICAL_COOLING/BG_TEMP_400K/FFSCT_0.100",
        # ROOT_PATH + "0.001_EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.100",
        # ROOT_PATH + "0.001_EMPIRICAL_COOLING/BG_TEMP_80K/FFSCT_0.100",
        # ROOT_PATH + "0.001_EMPIRICAL_COOLING/BG_TEMP_400K/FFSCT_0.100",
        # ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_10K/FFSCT_0.100/m_1.3_n_2.3",
        # ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/FFSCT_0.100/m_1.3_n_2.3",
        # ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_400K/FFSCT_0.100/m_1.3_n_2.3",
        # ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/FFSCT_0.100/m_1.5_n_3.0",
        # ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/FFSCT_0.047/m_1.3_n_2.3",
        # ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/FFSCT_0.100/m_2.0_n_4.0",
    # ]
    #
    # for target_folder_path in TARGET_FOLDER_PATHS:
    #     if not os.path.exists(target_folder_path):
    #         raise Exception("Folder does not exist!")
    #
    # timestamps = [
    #     # [1, 2, 7],
    #     # [1, 4, 11],
    #     # [1, 6, 14],
    #     [1, 8, 20],
    #     # [1, 4, 7],
    #     # [1, 6, 9],
    #     # [1, 2, 4],
    #     # [1, 4, 7],
    #     # [1, 5, 9],
    #     # [1, 6, 12],
    #     # [1, 6, 12],
    #     # [1, 6, 12],
    #     # [1, 3, 10],
    #     # [1, 3, 9],
    #     # [1, 3, 8],
    #     # [1, 2, 10],
    #     # [1, 2, 6],
    #     # [1, 2, 3],
    #     # [1, 2, 3],
    #     # [1, 2, 3],
    #     # [1, 2, 4],
    #     # [1, 2, 5],
    #     # [1, 2, 11]
    # ]

    main(readdir, writedir, timestamps)
    #To run compare diff eta


