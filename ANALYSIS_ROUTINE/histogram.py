import matplotlib
matplotlib.use('Agg')

import os
import argparse
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LogNorm, SymLogNorm, Normalize
import matplotlib.pyplot as plt
from pymses.sources.ramses.filename_utils import search_valid_outputs
from pymses.utils import constants as C
from pymses.filters import CellsToPoints
from pymses.analysis import Camera, slicing, ScalarOperator
import pymses

G0 = 1./1.7

def local_dir(read_directory, write_directory):
    if os.path.exists(read_directory):
        if not os.path.exists(write_directory):
            os.makedirs(write_directory)
    return read_directory, write_directory


def set_units(data):

    print('====== Output Info ======')
    for key in data.info:
        if key == "dom_decomp_Hilbert_keys":
            continue
        print(key, data.info[key])
    print ""

    boxlen = data.info["boxlen"]
    mu = 2.31  # 2.31  #ramses2.31, pymses 1/0.76
    fac_mu = 2.31 * 0.76
    fac_len = boxlen/data.info["unit_length"].express(C.pc)
    data.info["unit_length"] = data.info["unit_length"]*fac_len
    data.info["unit_density"] = data.info["unit_density"]
    data.info["unit_mass"] = data.info['unit_density'] * \
        data.info['unit_length'] ** 3
    #data.info["unit_mass"] = data.info["unit_mass"]*fac_len**3
    data.info["unit_time"] = data.info["unit_time"]
    data.info["unit_velocity"] = data.info["unit_velocity"]*fac_len
    data.info["unit_pressure"] = data.info["unit_pressure"]*fac_len**2
    data.info["unit_temperature"] = data.info["unit_temperature"] * \
        fac_len**2 * mu

    units['unit_M_g'] = data.info['unit_mass'].express(C.g)
    units['unit_M_kg'] = data.info['unit_mass'].express(C.kg)
    units['unit_M_Msun'] = data.info['unit_mass'].express(C.Msun)

    units['unit_L_cm'] = data.info['unit_length'].express(C.cm)
    units['unit_L_m'] = data.info['unit_length'].express(C.m)
    units['unit_L_km'] = data.info['unit_length'].express(C.km)
    units['unit_L_au'] = data.info['unit_length'].express(C.au)
    units['unit_L_pc'] = data.info['unit_length'].express(C.pc)

    units['unit_D_cgs'] = data.info['unit_density'].express(C.g/C.cm**3)
    units['unit_D_si'] = data.info['unit_density'].express(C.kg/C.m**3)
    units['unit_D_Hcc'] = data.info['unit_density'].express(
        C.H_cc) / 0.76 / 2.31
    #ts['unit_D_Hcc'] = data.info['unit_density'].express(C.H_cc)

    units['unit_T_s'] = data.info['unit_time'].express(C.s)
    units['unit_T_yr'] = data.info['unit_time'].express(C.year)
    units['unit_T_Kyr'] = data.info['unit_time'].express(C.year) * (10 ** -3)
    units['unit_T_Myr'] = data.info['unit_time'].express(C.year) * (10 ** -6)

    units['unit_V_cgs'] = data.info['unit_velocity'].express(C.cm/C.s)
    units['unit_V_si'] = data.info['unit_velocity'].express(C.m/C.s)

    units['unit_Pressure'] = data.info['unit_pressure'].express(C.N/C.m**2)
    units['unit_Temperature'] = data.info['unit_temperature'].express(C.K)

    print("====== Unit after normalize ======")
    for item in units.items():
        print(item)
    print ""
    return units


def plot_density_profile(ax, x, y):
    ax.plot(x, y, "-", markersize=2, color="black")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10**-2, 10**2)
    ax.set_ylim(10**0, 10**7)
    ax.set_xlabel("Radius (pc)")
    ax.set_ylabel("Number Density " + r"$\rm(1/cm^3)$")
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    return


def plot_temperature_profile(ax, x, y):
    ax.plot(x, y, "-", markersize=2, color="black")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10**-2, 10**2)
    ax.set_ylim(10**0, 10**6)
    ax.set_xlabel("Radius (pc)")
    ax.set_ylabel("Temperature " + r"$\rm(K)$")
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    return


def plot_radial_velocity_profile(ax, x, y):
    ax.plot(x, y, "-", markersize=2, color="black")
    ax.set_xscale("log")
    ax.set_xlim(10**-2, 10**2)
    ax.set_ylim(np.nanmin(y), np.nanmax(y))
    ax.set_xlabel("Radius (pc)")
    ax.set_ylabel("Radial Velocity " + r"$\rm(m/s)$")
    ax.axhline(y=0.0, c="r", ls="--", lw=1)
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    return


def plot_enclosed_mass(ax, x, y):
    ax.plot(x, y, "-", markersize=2, color="black")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10 ** -2, 10 ** 2)
    ax.set_xlabel("Radius (pc)")
    ax.set_ylabel("Enclosed Mass " + r"$\rm(M_{\odot})$")
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    return


def plot_density_temperature_2d_histogram_pcolor(fig, ax, hist_2d, density_edges, temperature_edges):
    im = ax.pcolormesh(density_edges, temperature_edges, hist_2d,
                       norm=LogNorm(), vmin=1e0, vmax=1e1, cmap='rainbow')
    ax.set_xlabel("Number Density " + r"$\rm(1/cm^3)$")
    ax.set_ylabel("Temperature " + r"$\rm(K)$")
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    fig.colorbar(im, ax=ax, extend='max', label='Mass ' + r'$\rm(M_{\odot})$')
    return


def plot_density_temperature_2d_histogram_imshow(fig, ax, hist_2d, density_edges, temperature_edges):

    im = ax.imshow(hist_2d, origin='lower', norm=LogNorm(), cmap='rainbow', vmin=1e-10, vmax=1e4,
                   extent=(density_edges[0], density_edges[-1], temperature_edges[0], temperature_edges[-1]))

    ax.set_xlabel("Log10 Number Density " + r"$\rm(1/cm^3)$")
    ax.set_ylabel("Log10 Temperature " + r"$\rm(K)$")
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    fig.colorbar(im, ax=ax, extend='both', label='Mass ' + r'$\rm(M_{\odot})$')
    return


def plot_theoretical_equilibrium_temperature(ax, x, y):
    ax.plot(x, y, '-', lw=2, color='black', alpha=0.2)
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    return


def make_histogram(read_dir, write_dir, output):

    global units
    bins = np.logspace(-4.0, -0.3, 32)
    bins_center = np.sqrt(bins[:-1] * bins[1:])

    bins = np.append(0, bins)
    bins_center = np.append(bins[1]/2, bins_center)

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

    # Calculate density profile
    mass_histogram, _ = np.histogram(radius, bins=bins, weights=mass)
    volume_histogram, _ = np.histogram(radius, bins=bins, weights=volume)
    density_histogram = mass_histogram / volume_histogram

    # Calculate temperature profile
    mass_weighted_temperature_histogram, _ = np.histogram(
        radius, bins=bins, weights=temperature * mass)
    temperature_histogram = mass_weighted_temperature_histogram / mass_histogram

    # Calculate radial velocity profile
    mass_weighted_radial_velocity_histogram, _ = np.histogram(
        radius, bins=bins, weights=radial_velocity * mass)
    radial_velocity_histogram = mass_weighted_radial_velocity_histogram / mass_histogram

    # Calculate enclosed mass profile
    enclosed_mass = np.cumsum(mass_histogram)

    # Calculate Density vs. Temperature 2D histogram
    density_bins = np.linspace(-2, 8, 128)
    temparture_bins = np.linspace(-2, 8, 128)
    density_temperature_2d_histogram, density_edges, temperature_edges = np.histogram2d(
        np.log10(density), np.log10(temperature), weights=mass, bins=(density_bins, temparture_bins))
    density_temperature_2d_histogram = np.where(
        density_temperature_2d_histogram > 0, density_temperature_2d_histogram, np.nan)

    # Calculate theoretical equilibrium temperature
    nH = np.load('/data/daniellin/COOLING_FUNCTION_ANALYZE/EQUILIBRIUM/nH.npy')
    equilibrium_temperature = np.load(
        '/data/daniellin/COOLING_FUNCTION_ANALYZE/EQUILIBRIUM/equilibrium_Temp.npy')

    # Plot Figures
    fig, ((ax0, ax1, ax2), (ax3, ax4, ax5)) = plt.subplots(
        nrows=2, ncols=3, figsize=(21, 12))
    fig.suptitle("No. " + str(output) + " Time: " +
                 str(round(data.info["time"] * units['unit_T_Kyr'], 3)) + " (kyr)")

    plot_density_profile(
        ax0, bins_center * units['unit_L_pc'], density_histogram * units["unit_D_Hcc"])
    plot_temperature_profile(
        ax1, bins_center * units['unit_L_pc'], temperature_histogram)
    plot_radial_velocity_profile(
        ax2, bins_center * units['unit_L_pc'], radial_velocity_histogram * units["unit_V_si"])
    plot_enclosed_mass(
        ax3, bins_center * units['unit_L_pc'], enclosed_mass * units["unit_M_Msun"])

    
    plot_density_temperature_2d_histogram_imshow(fig, ax4, density_temperature_2d_histogram.T *
        units["unit_M_Msun"], density_edges + np.log10(units["unit_D_Hcc"]), temperature_edges)
    plot_theoretical_equilibrium_temperature(
        ax4, np.log10(nH), np.log10(equilibrium_temperature))

    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.99], pad=1.8)
    plt.savefig(write_dir + '/profiles_' +
                format(output, '05') + ".pdf", bbox_inches='tight')
    #plt.show()
    plt.clf()
    return


units = dict()
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="plot profiles from RAMSES output")
    parser.add_argument('-i', '--input', type=str, help="RAMSES output path")
    parser.add_argument('-o', '--output', type=str, help="Figure output path")
    args = parser.parse_args()

    read_directory, write_directory = local_dir(args.input, args.output)
    outputs = search_valid_outputs(read_directory)

    for output in outputs[3:]:
        make_histogram(read_directory, write_directory, output)
