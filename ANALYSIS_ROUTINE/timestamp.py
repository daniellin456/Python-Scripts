import os
import pymses
import csv
import pandas as pd
import numpy as np
from pymses.sources.ramses.filename_utils import  search_valid_outputs
from module_units import *

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

def get_parametric_cooling_m_n(target_folder_path):
    m = ""
    n = ""
    if "PARAMETRIC_COOLING" in target_folder_path:
        if "m_" in target_folder_path:
            index = target_folder_path.index("m_")
            m = target_folder_path[index+2:index+5]
            index = target_folder_path.index("n_")
            n = target_folder_path[index+2:index+5]
    else:
        return

    return m, n

if __name__ == "__main__":
    ROOT_PATH = "/data/daniellin/RAMSES_RUNS/"
    TARGET_FOLDER_PATHS = [
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
        # ROOT_PATH + "0.001_EMPIRICAL_COOLING/BG_TEMP_400K/FFSCT_0.100"
        ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_10K/FFSCT_0.100/m_1.3_n_2.3",
        ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/FFSCT_0.100/m_1.3_n_2.3",
        ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_400K/FFSCT_0.100/m_1.3_n_2.3",
        ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/FFSCT_0.100/m_1.5_n_3.0",
        ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/FFSCT_0.047/m_1.3_n_2.3",
        ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/FFSCT_0.100/m_2.0_n_4.0"
    ]

    df = pd.DataFrame(columns=["timestamp"])

    for target_folder_path in TARGET_FOLDER_PATHS:
        if not os.path.exists(target_folder_path):
            raise Exception("Folder does not exist!")

        print(target_folder_path)
        outputs = search_valid_outputs(target_folder_path)
        eta = find_eta(target_folder_path)
        T0 = find_T0(target_folder_path)
        ffsct = find_ffsct(target_folder_path)
        m, n = get_parametric_cooling_m_n(target_folder_path)

        if "PARAMETRIC_COOLING" in target_folder_path:
            title = "m_" + m + "_n_" + n + "_T0_" + T0 + "_ffsct_" + ffsct
            print(m, n)
        else:
            title = "eta_" + eta + "_T0_" + T0 + "_ffsct_" + ffsct
        df[title] = ""

        for output in outputs[0:len(outputs)-1]:
            data = pymses.RamsesOutput(target_folder_path, output)
            free_fall_time = find_free_fall_time(target_folder_path)
            units = set_units(data)
            time_sec = data.info['time'] * units['unit_T_s']
            time_kyr = data.info['time'] * units['unit_T_Kyr']
            ratio_to_free_fall_time = round(time_sec / free_fall_time, 3)
            df.loc[output, title] = ratio_to_free_fall_time

    df.to_csv("output_parametric.csv", index=False)