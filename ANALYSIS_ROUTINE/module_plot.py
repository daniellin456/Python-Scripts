import numpy as np
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LogNorm, SymLogNorm, Normalize
import matplotlib.pyplot as plt
import numpy as np


def plot_density_profile(ax, x, y):
    ax.plot(x, y, "-", markersize=2, color="black")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10 ** -2, 10 ** 2)
    ax.set_ylim(10 ** 0, 10 ** 7)
    ax.set_xlabel("Radius (pc)")
    ax.set_ylabel("Density " + r"$\rm(cm^{-3})$")
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    return


def plot_temperature_profile(ax, x, y):
    ax.plot(x, y, "-", markersize=2, color="black")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10 ** -2, 10 ** 2)
    ax.set_ylim(10 ** 0, 10 ** 6)
    ax.set_xlabel("Radius (pc)")
    ax.set_ylabel("Temperature " + r"$\rm(K)$")
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    return


def plot_radial_velocity_profile(ax, x, y):
    ax.plot(x, y, "-", markersize=2, color="black")
    ax.set_xscale("log")
    ax.set_xlim(10 ** -2, 10 ** 2)
    ax.set_ylim(np.nanmin(y), np.nanmax(y))
    ax.set_xlabel("Radius (pc)")
    ax.set_ylabel("Radial Velocity " + r"$\rm(ms^{-1})$")
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

    # ax.set_xlabel(r"$Log10 Density " + r"$\rm(1/cm^3)$")
    ax.set_xlabel(r"$\rm{log_{10}}\; n_H \; (\rm{cm^{-3}})$")
    # ax.set_ylabel(r"$Log10 Temperature " + r"$\rm(K)$")
    ax.set_ylabel(r"$\rm{log_{10}}\; T \; (\rm{K})$")
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    fig.colorbar(im, ax=ax, extend='both', orientation="horizontal", label='Mass ' + r'$\rm(M_{\odot})$')
    return


def plot_theoretical_equilibrium_temperature(ax, x, y):
    ax.plot(x, y, '--', lw=2, color='black')
    ax.tick_params(which='major', width=1, length=5)
    ax.tick_params(which='minor', width=1, length=3)
    return


def plot_slope_line(ax, slope, intercepts, color):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_vals = np.array(xlim)
    for intercept in intercepts:
        y_vals = intercept + slope * x_vals
        ax.plot(x_vals, y_vals, linestyle=':', color=color)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return


def plot_integrated(data, write_dir, output, units, nH, bins_center, density_edges, temperature_edges,
                    density_histogram,
                    temperature_histogram, radial_velocity_histogram, enclosed_mass, density_temperature_2d_histogram,
                    equilibrium_temperature):
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
                                                 units["unit_M_Msun"], density_edges + np.log10(units["unit_D_Hcc"]),
                                                 temperature_edges)

    plot_slope_line(ax4, 0.666667, np.arange(-6, 4, 2), 'black')

    plot_theoretical_equilibrium_temperature(
        ax4, np.log10(nH), np.log10(equilibrium_temperature))

    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.99], pad=1.8)
    plt.savefig(write_dir + '/profiles_' +
                format(output, '05') + ".pdf", bbox_inches='tight')
    return


def plot_single(data, write_dir, output, units, nH, bins_center, density_edges, temperature_edges, density_histogram,
                temperature_histogram, radial_velocity_histogram, enclosed_mass, density_temperature_2d_histogram,
                equilibrium_temperature):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 6))
    ax.set_title("Time: " + str(round(data.info["time"] * units['unit_T_Kyr'], 3)) + " (kyr)")
    # ax.set_title("Density histogram")
    plot_density_profile(ax, bins_center * units['unit_L_pc'], density_histogram * units["unit_D_Hcc"])
    plt.savefig(write_dir + '/dens_' + format(output, '05') + ".pdf", bbox_inches='tight')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 6))
    ax.set_title("Time: " + str(round(data.info["time"] * units['unit_T_Kyr'], 3)) + " (kyr)")
    # ax.set_title("Temperature histogram")
    plot_temperature_profile(ax, bins_center * units['unit_L_pc'], temperature_histogram)
    plt.savefig(write_dir + '/temp_' + format(output, '05') + ".pdf", bbox_inches='tight')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 6))
    ax.set_title("Time: " + str(round(data.info["time"] * units['unit_T_Kyr'], 3)) + " (kyr)")
    # ax.set_title("Radial velocity histogram")
    plot_radial_velocity_profile(ax, bins_center * units['unit_L_pc'], radial_velocity_histogram * units["unit_V_si"])
    plt.savefig(write_dir + '/vel_' + format(output, '05') + ".pdf", bbox_inches='tight')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 6))
    ax.set_title("Time: " + str(round(data.info["time"] * units['unit_T_Kyr'], 3)) + " (kyr)")
    # ax.set_title("Enclosed mass")
    plot_enclosed_mass(ax, bins_center * units['unit_L_pc'], enclosed_mass * units["unit_M_Msun"])
    plt.savefig(write_dir + '/mass_' + format(output, '05') + ".pdf", bbox_inches='tight')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 6))
    ax.set_title("Time: " + str(round(data.info["time"] * units['unit_T_Kyr'], 3)) + " (kyr)")
    # ax.set_title("Density-Temperature relation")
    plot_density_temperature_2d_histogram_imshow(fig, ax, density_temperature_2d_histogram.T *
                                                 units["unit_M_Msun"], density_edges + np.log10(units["unit_D_Hcc"]),
                                                 temperature_edges)
    plot_theoretical_equilibrium_temperature(ax, np.log10(nH), np.log10(equilibrium_temperature))
    plt.savefig(write_dir + '/dens_temp_' + format(output, '05') + ".pdf", bbox_inches='tight')
    return
