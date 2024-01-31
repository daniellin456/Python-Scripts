import os
import pymses
import csv
import pandas as pd
import numpy as np
from pymses.sources.ramses.filename_utils import  search_valid_outputs
from module_units import *

def find_free_fall_time(target_folder):

    peak_density = 0.0
    if "SC01" in target_folder:
        peak_density = 1280.01295051676

    elif "SC0068" in target_folder:
        peak_density = 12946.7537409113

    elif "SC0047" in target_folder:
        peak_density = 118748.249427123

    free_fall_time = np.sqrt(3.0 * np.pi / 32.0 / 6.67e-8 / 2.31 / 1.66e-24 / peak_density)

    return free_fall_time

if __name__ == "__main__":
    max_output = 0
    ROOT_PATH = "/data/daniellin/RAMSES_RUNS/"
    folders = os.listdir(ROOT_PATH)
    df = pd.DataFrame(columns=["timestamp"])

    for folder in folders:
        target_folder = ROOT_PATH + folder
        print(target_folder)
        outputs = search_valid_outputs(target_folder)
        df[folder] = ""

        for output in outputs[0:len(outputs)-1]:
            data = pymses.RamsesOutput(target_folder, output)
            free_fall_time = find_free_fall_time(target_folder)
            units = set_units(data)
            time_sec = data.info['time'] * units['unit_T_s']
            time_kyr = data.info['time'] * units['unit_T_Kyr']
            ratio_to_free_fall_time = round(time_sec / free_fall_time, 3)
            df.loc[output, folder] = ratio_to_free_fall_time

        for i in range(1, len(outputs)):
            df.loc[i, "timestamp"] = i

    df.to_csv("output_parametric.csv", index=False)
