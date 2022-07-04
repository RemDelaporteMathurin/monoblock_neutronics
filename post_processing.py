import openmc
import regular_mesh_plotter as rmp
import matplotlib.pyplot as plt
import openmc_tally_unit_converter as otuc
from scipy.stats import linregress
import numpy as np


def shape_tally(tally, threshold=0):
    value = np.array(tally.mean)
    std_dev = np.array(tally.std_dev)

    voxel_volume = otuc.compute_volume_of_voxels(tally)  # cm3
    voxel_volume *= 1e-6  # m3
    value *= 1 / voxel_volume  # m-3 per neutron
    value *= source_strength  # m-3 s-1

    std_dev *= 1 / voxel_volume  # m-3 per neutron
    std_dev *= source_strength  # m-3 s-1

    value = rmp.reshape_values_to_mesh_shape(tally, value)
    std_dev = rmp.reshape_values_to_mesh_shape(tally, std_dev)

    indexes = np.where(value < threshold)

    value[indexes] = np.nan
    std_dev[indexes] = np.nan
    return value, std_dev


def plot_distribution(mean, std_dev):
    distribution = np.nanmean(mean, axis=1)
    std_dev_distrib = np.nanmean(std_dev, axis=1)
    std_dev_distrib = std_dev_distrib[np.logical_not(np.isnan(distribution))]
    distribution = distribution[np.logical_not(np.isnan(distribution))]
    x = np.linspace(0, height_mb, len(distribution), endpoint=True)

    res = linregress(x[1:], distribution[::-1][1:])

    plt.errorbar(
        x, distribution[::-1], yerr=2 * std_dev_distrib[::-1], fmt="o", alpha=0.5
    )
    (line,) = plt.plot(x, res.slope * x + res.intercept)
    x_annotation = 0.75
    y_annotation = (res.slope * x_annotation + res.intercept) * 1.1
    plt.annotate(
        "{:.1e} $x$ + {:.1e}".format(res.slope, res.intercept),
        (x_annotation, y_annotation),
        color=line.get_color(),
    )


result = openmc.StatePoint("statepoint.50.h5")

height_mb = 2.5  # cm
source_strength = otuc.find_source_strength(
    fusion_energy_per_second_or_per_pulse=500e6
)  # n/s

# He generation mesh tally
helium_generation_mesh = result.get_tally(name="(n,Xa)_on_2D_mesh_yz")

value, std_dev = shape_tally(helium_generation_mesh, threshold=4.75e18)

rmp.plot_regular_mesh_values(
    values=value,
    extent=rmp.get_tally_extent(helium_generation_mesh),
    rotate_plot=180,
    label="He3 generation (m$^{-3}$ s$^{-1}$)",
)
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.savefig("helium_generation.png")

# depth distribution of helium generation

plt.figure()
plot_distribution(value, std_dev)
plt.xlabel("Distance from the top surface (cm)")
plt.ylabel("He generation  (m$^{-3}$ s$^{-1}$)")
plt.ylim(bottom=0)
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.savefig("he_generation_distribution.png")


# Heat generation mesh tally
plt.figure()
heat_generation_mesh = result.get_tally(name="heating_on_2D_mesh_yz")

value, std_dev = shape_tally(heat_generation_mesh, threshold=4.22e27)
value *= 1.602e-19 * 1e-6  # MW  m-3
std_dev *= 1.602e-19 * 1e-6  # MW  m-3

rmp.plot_regular_mesh_values(
    values=value,
    extent=rmp.get_tally_extent(heat_generation_mesh),
    rotate_plot=180,
    label="Heat generation (MW m$^{-3}$)",
)
plt.gca().get_images()[0].set_cmap("inferno")
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.savefig("heat_generation.png")

# depth distribution of helium generation
plt.figure()

plot_distribution(value, std_dev)
plt.xlabel("Distance from the top surface (cm)")
plt.ylabel("Heat generation  (MW m$^{-3}$)")
plt.ylim(bottom=0)
plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.savefig("heat_distribution.png")
