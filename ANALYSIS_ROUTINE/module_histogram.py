import numpy as np


def calculate_mass_histogram(radius, bins, weights):
    mass_histogram, _ = np.histogram(radius, bins=bins, weights=weights)
    return mass_histogram


def calculate_volume_histogram(radius, bins, weights):
    volume_histogram, _ = np.histogram(radius, bins=bins, weights=weights)
    return volume_histogram


def calculate_density_histogram(mass_histogram, volume_histogram):
    return mass_histogram / volume_histogram


def calculate_temperature_histogram(radius, bins, weights, mass_histogram):
    mass_weighted_temperature_histogram, _ = np.histogram(
        radius, bins=bins, weights=weights[0] * weights[1])
    temperature_histogram = mass_weighted_temperature_histogram / mass_histogram
    return temperature_histogram


def calculate_radial_velocity_histogram(radius, bins, weights, mass_histogram):
    mass_weighted_radial_velocity_histogram, _ = np.histogram(
        radius, bins=bins, weights=weights[0] * weights[1])
    radial_velocity_histogram = mass_weighted_radial_velocity_histogram / mass_histogram
    return radial_velocity_histogram


def calculate_enclosed_mass(mass_histogram):
    return np.cumsum(mass_histogram)


def calculate_density_temperature_2d_histogram(density, temperature, weights):
    density_bins = np.linspace(-2, 10, 128)
    temparture_bins = np.linspace(0, 4, 128)
    density_temperature_2d_histogram, density_edges, temperature_edges = np.histogram2d(
        np.log10(density), np.log10(temperature), weights=weights, bins=(density_bins, temparture_bins))
    density_temperature_2d_histogram = np.where(
        density_temperature_2d_histogram > 0, density_temperature_2d_histogram, np.nan)

    return density_temperature_2d_histogram, density_edges, temperature_edges


def calculate_equilibrium_temperature():
    nH = np.load('/data/daniellin/PYTHON_SCRIPTS/COOLING_FUNCTION_ANALYZE/EQUILIBRIUM/nH.npy')
    equilibrium_temperature = np.load(
        '/data/daniellin/PYTHON_SCRIPTS/COOLING_FUNCTION_ANALYZE/EQUILIBRIUM/equilibrium_Temp.npy')
    return nH, equilibrium_temperature
