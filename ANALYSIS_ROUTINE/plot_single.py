import os
import matplotlib

matplotlib.use('Agg')

import argparse
from pymses.sources.ramses.filename_utils import search_valid_outputs
from pymses.filters import CellsToPoints
import pymses
from module_plot import *
from module_histogram import *
from module_units import *


def local_dir(read_directory, write_directory):
    if os.path.exists(read_directory):
        if not os.path.exists(write_directory):
            os.makedirs(write_directory)
    return read_directory, write_directory


def make_histogram(read_dir, write_dir, output):
    global units

    bins = np.logspace(-4.0, -0.3, 32)
    bins_center = np.sqrt(bins[:-1] * bins[1:])

    bins = np.append(0, bins)
    bins_center = np.append(bins[1] / 2, bins_center)

    data = pymses.RamsesOutput(read_dir, output)

    if len(units) == 0:
        units = set_units(data)

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
    density_temperature_2d_histogram, density_edges, temperature_edges = calculate_density_temperature_2d_histogram(
        density, temperature, mass)
    nH, equilibrium_temperature = calculate_equilibrium_temperature()
    plot_single(data, write_dir, output, units, nH, bins_center, density_edges, temperature_edges, density_histogram,
                temperature_histogram, radial_velocity_histogram, enclosed_mass, density_temperature_2d_histogram,
                equilibrium_temperature)
    return


units = dict()
G0 = 1. / 1.7
if __name__ == '__main__':
    ROOT_PATH = "/data/daniellin/RAMSES_RUNS/"
    target_folders = [#ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.100",
                      #ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.068",
                      #ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_10K/FFSCT_0.047",
                      #ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_80K/FFSCT_0.100",
                      #ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_80K/FFSCT_0.068",
                      #ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_80K/FFSCT_0.047",
                      #ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_400K/FFSCT_0.100",
                      #ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_400K/FFSCT_0.068",
                      #ROOT_PATH + "EMPIRICAL_COOLING/BG_TEMP_400K/FFSCT_0.047",
                      #ROOT_PATH + "0.1_EMPIRICAL_COOLING/BG_TEMP_10K",
                      #ROOT_PATH + "0.1_EMPIRICAL_COOLING/BG_TEMP_80K",
                      #ROOT_PATH + "0.1_EMPIRICAL_COOLING/BG_TEMP_400K",
                      #ROOT_PATH + "0.01_EMPIRICAL_COOLING/BG_TEMP_10K",
                      #ROOT_PATH + "0.01_EMPIRICAL_COOLING/BG_TEMP_80K",
                      #ROOT_PATH + "0.01_EMPIRICAL_COOLING/BG_TEMP_400K",
                      #ROOT_PATH + "0.001_EMPIRICAL_COOLING/BG_TEMP_10K",
                      #ROOT_PATH + "0.001_EMPIRICAL_COOLING/BG_TEMP_80K",
                      #ROOT_PATH + "0.001_EMPIRICAL_COOLING/BG_TEMP_400K",
                      ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/m_1.3_n_2.3",
                      ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/m_1.5_n_3.0",
                      ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_80K/m_2.0_n_4.0",
                      ROOT_PATH + "PARAMETRIC_COOLING/BG_TEMP_400K",
                      ROOT_PATH + "PARAMETRIC_COOLING/FFSCT_0.047/m_1.3_n_2.3"]

    # parser = argparse.ArgumentParser(
    #     description="plot profiles from RAMSES output")
    # parser.add_argument('-i', '--input', type=str, help="RAMSES output path")
    # parser.add_argument('-o', '--output', type=str, help="Figure output path")
    # args = parser.parse_args()
    #
    # read_directory, write_directory = local_dir(args.input, args.output)

    for folders in target_folders:
        units.clear()
        outputs = search_valid_outputs(folders)

        for output in outputs[:3]:
            make_histogram(folders, folders+"/profiles", output)
