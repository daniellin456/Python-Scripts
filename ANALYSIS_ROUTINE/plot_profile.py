import os
import argparse
import matplotlib
matplotlib.use('Agg')

import pymses
from pymses.sources.ramses.filename_utils import  search_valid_outputs
from pymses.filters import CellsToPoints
from module_style import *
from module_histogram import *
from module_units import *
from module_plot import *

def find_free_fall_time(target_folder_path):

    peak_density = 0.0
    if "0.100" in target_folder_path:
        peak_density = 1280.01295051676

    elif "0.068" in target_folder_path:
        peak_density = 12946.7537409113

    elif "0.047" in target_folder_path:
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
    if "0.100" in target_folder_path:
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

def get_density_xlim_ylim(ffsct):
    xlim = []
    ylim = []
    if ffsct == "0.100":
        xlim = [10**-1.4, 10**1.3]
        ylim = [10**0.5, 10**4.0]
    elif ffsct == "0.068":
        xlim = [10 ** -1.4, 10 ** 1.0]
        ylim = [10 ** 1.5, 10 ** 5.0]
    elif ffsct == "0.047":
        xlim = [10 ** -1.7, 10 ** 0.5]
        ylim = [10 ** 2.5, 10 ** 5.5]
    return xlim, ylim

def get_temperature_xlim_ylim(T0):
    xlim = []
    ylim = []
    if T0 == "10K":
        ylim = [10 ** 0.5, 10 ** 3.0]
    elif T0 == "80K":
        ylim = [10 ** 0.5, 10 ** 3.0]
    elif T0 == "400K":
        ylim = [10 ** 0.5, 10 ** 3.0]
    return ylim

def get_parametric_cooling_m_n(target_folder_path):
    m = ""
    n = ""
    if "PARAMETRIC_COOLING" in target_folder_path:
        if "m_" in target_folder_path:
            index = target_folder_path.index("m_")
            m = target_folder_path[index + 2:index + 5]
            index = target_folder_path.index("n_")
            n = target_folder_path[index + 2:index + 5]
    else:
        return m, n

    return m, n

def main(target_folder_paths, timestamps):

    for target_folder_path in target_folder_paths: # Loop FFSCT = 0.1, 0.068, 0.047
        fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(7.0, 7.0))
        fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(7.0, 7.0))
        fig3, ax3 = plt.subplots(nrows=1, ncols=1, figsize=(7.0, 7.0))
        fig4, ax4 = plt.subplots(nrows=1, ncols=1, figsize=(7.0, 7.0))
        free_fall_time = find_free_fall_time(target_folder_path)
        ffsct = find_ffsct(target_folder_path)
        T0 = find_T0(target_folder_path)
        eta = find_eta(target_folder_path)
        m, n = get_parametric_cooling_m_n(target_folder_path)
        color_list = ['red', 'green', 'blue']

        for timestamp, color in zip(timestamps[target_folder_paths.index(target_folder_path)], color_list): # Draw 3 different timestamp in each FFSCT

            bins = np.logspace(-2.0, -0.3, 32)
            bins_center = np.sqrt(bins[:-1] * bins[1:])
            bins = np.append(0, bins)
            bins_center = np.append(bins[1] / 2, bins_center)

            data = pymses.RamsesOutput(target_folder_path, timestamp)
            units = set_units(data)
            time_sec = data.info['time'] * units['unit_T_s']
            time_kyr = data.info['time'] * units['unit_T_Kyr']
            amr_data = data.amr_source(["rho", "T", "vel"])
            cell = CellsToPoints(amr_data).flatten()
            coordinate = cell.points
            center_coordinate = np.array([0.5, 0.5, 0.5])
            volume = cell.get_sizes() ** 3
            density = cell["rho"]
            temperature = cell["T"]
            velocity = cell["vel"]
            mass = density * volume
            center_coordinate = center_coordinate[np.newaxis, :]
            radius = np.linalg.norm(coordinate - center_coordinate, axis=1)
            radial_velocity = np.sum(np.multiply(
                velocity, coordinate - center_coordinate), axis=1) / radius

            mass_histogram = calculate_mass_histogram(radius, bins, mass)
            volume_histogram = calculate_volume_histogram(radius, bins, volume)
            density_histogram = calculate_density_histogram(mass_histogram, volume_histogram)
            temperature_histogram = calculate_temperature_histogram(radius, bins, [temperature, mass], mass_histogram)
            radial_velocity_histogram = calculate_radial_velocity_histogram(radius, bins, [radial_velocity, mass],
                                                                            mass_histogram)
            enclosed_mass = calculate_enclosed_mass(mass_histogram)
            legend = "t = " + str(round(time_sec/free_fall_time, 3)) + r"$ t_{ff}$, " + str(round(time_kyr, 3)) + " kyr"


            plot_density_profile(
                ax1, bins_center * units['unit_L_pc'], density_histogram * units["unit_D_Hcc"], color=color, label=legend)
            plot_temperature_profile(
                ax2, bins_center * units['unit_L_pc'], temperature_histogram, color=color, label=legend)
            plot_radial_velocity_profile(
                ax3, bins_center * units['unit_L_pc'], radial_velocity_histogram * units["unit_V_si"], color=color, label=legend)
            plot_enclosed_mass(
                ax4, bins_center * units['unit_L_pc'], enclosed_mass * units["unit_M_Msun"], color=color, label=legend)

        if m != "" and n != "":
            ax1.set_title(r"$t_{ff}/t_{sc}$ = " + ffsct + ", m = " + m + ", n = " + n+ r", $T_0$ = " + T0)#, fontsize = 8)
            ax2.set_title(r"$t_{ff}/t_{sc}$ = " + ffsct + ", m = " + m + ", n = " + n+ r", $T_0$ = " + T0)#, fontsize = 8)
            ax3.set_title(r"$t_{ff}/t_{sc}$ = " + ffsct + ", m = " + m + ", n = " + n+ r", $T_0$ = " + T0)#, fontsize = 8)
            ax4.set_title(r"$t_{ff}/t_{sc}$ = " + ffsct + ", m = " + m + ", n = " + n+ r", $T_0$ = " + T0)#, fontsize = 8)
        else:
            ax1.set_title(r"$t_{ff}/t_{sc}$ = " + ffsct + r", $\eta$ = " + eta + r", $T_0$ = " + T0)
            ax2.set_title(r"$t_{ff}/t_{sc}$ = " + ffsct + r", $\eta$ = " + eta + r", $T_0$ = " + T0)
            ax3.set_title(r"$t_{ff}/t_{sc}$ = " + ffsct + r", $\eta$ = " + eta + r", $T_0$ = " + T0)
            ax4.set_title(r"$t_{ff}/t_{sc}$ = " + ffsct + r", $\eta$ = " + eta + r", $T_0$ = " + T0)

        ax1.legend(loc='upper right')#, fontsize=7)
        ax2.legend(loc='upper right')#, fontsize=7)
        ax3.legend(loc='upper right')#, fontsize=7)
        ax4.legend(loc='upper right')#, fontsize=7)
        fig1.tight_layout()
        fig2.tight_layout()
        fig3.tight_layout()
        fig4.tight_layout()

        if m != "" and n != "":
            fig1.savefig("m_" + m + "_n_" + n + "_T0_" + T0 + "_ffsct_" + ffsct + "_Dens.pdf")
            fig2.savefig("m_" + m + "_n_" + n + "_T0_" + T0 + "_ffsct_" + ffsct + "_Temp.pdf")
            fig3.savefig("m_" + m + "_n_" + n + "_T0_" + T0 + "_ffsct_" + ffsct + "_Vel.pdf")
            fig4.savefig("m_" + m + "_n_" + n + "_T0_" + T0 + "_ffsct_" + ffsct + "_Mass.pdf")
        else:
            fig1.savefig("eta_" + eta + "_T0_" + T0 + "_ffsct_" + ffsct + "_Dens.pdf")
            fig2.savefig("eta_" + eta + "_T0_" + T0 + "_ffsct_" + ffsct + "_Temp.pdf")
            fig3.savefig("eta_" + eta + "_T0_" + T0 + "_ffsct_" + ffsct + "_Vel.pdf")
            fig4.savefig("eta_" + eta + "_T0_" + T0 + "_ffsct_" + ffsct + "_Mass.pdf")
    return

if __name__ == "__main__":

    plot_style()

    ROOT_PATH = "/data/daniellin/RAMSES_RUNS/"
    TARGET_FOLDER_PATHS = [
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.100",
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.068",
        # ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.047",
        ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_80K/FFSCT_0.100",
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
        # ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/FFSCT_0.100/m_2.0_n_4.0"
    ]

    for target_folder_path in TARGET_FOLDER_PATHS:
        if not os.path.exists(target_folder_path):
            raise Exception("Folder does not exist!")

    timestamps = [
        # [1, 2, 7],
        # [1, 4, 11],
        # [1, 6, 14],
        [1, 8, 20],
        # [1, 4, 7],
        # [1, 6, 9],
        # [1, 2, 4],
        # [1, 4, 7],
        # [1, 5, 9],
        # [1, 6, 12],
        # [1, 6, 12],
        # [1, 6, 12],
        # [1, 3, 10],
        # [1, 3, 9],
        # [1, 3, 8],
        # [1, 2, 10],
        # [1, 2, 6],
        # [1, 2, 3],
        # [1, 2, 3],
        # [1, 2, 3],
        # [1, 2, 4],
        # [1, 2, 5],
        # [1, 2, 11]
    ]

    main(TARGET_FOLDER_PATHS, timestamps)
    #To run compare diff eta


