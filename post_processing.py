import openmc
import regular_mesh_plotter as rmp
import matplotlib.pyplot as plt
import openmc_tally_unit_converter as otuc
from scipy.stats import linregress
import numpy as np

result = openmc.StatePoint("statepoint.20.h5")

source_strength = otuc.find_source_strength(
    fusion_energy_per_second_or_per_pulse=500e6
)  # n/s

total_helium_generation_in_tungsten = result.get_tally(name="tungsten_(n,Xa)")

# He generation mesh tally
helium_generation_mesh = result.get_tally(name="(n,Xa)_on_2D_mesh_yz")

value = helium_generation_mesh.mean

voxel_volume = otuc.compute_volume_of_voxels(helium_generation_mesh)  # cm3
voxel_volume *= 1e-6  # m3


value *= 1 / voxel_volume  # m-3 per neutron
value *= source_strength  # m-3 s-1

value = rmp.reshape_values_to_mesh_shape(helium_generation_mesh, value)
threshold = 0.6e19
value[np.where(value < threshold)] = np.nan

rmp.plot_regular_mesh_values(
    values=value,
    extent=rmp.get_tally_extent(helium_generation_mesh),
    rotate_plot=180,
    label="He3 generation (m$^{-3}$ s$^{-1}$)",
)
plt.savefig("helium_generation.png")

# depth distribution of helium generation
distribution = np.nanmean(value, axis=1)
distribution = distribution[np.logical_not(np.isnan(distribution))]
height_mb = 2.5  # cm
x = np.linspace(0, height_mb, len(distribution), endpoint=True)

res = linregress(x[1:], distribution[::-1][1:])

plt.figure()
plt.scatter(x, distribution[::-1])
plt.plot(x, res.slope * x + res.intercept)
plt.xlabel("Distance from the top surface (cm)")
plt.ylabel("He generation  (m$^{-3}$ s$^{-1}$)")
plt.ylim(bottom=0, top=7e18)
plt.savefig("he_generation_distribution.png")
# Heat generation mesh tally
plt.figure()
heat_generation_mesh = result.get_tally(name="heating_on_2D_mesh_yz")

value = heat_generation_mesh.mean

voxel_volume = otuc.compute_volume_of_voxels(heat_generation_mesh)  # cm3
voxel_volume *= 1e-6  # m3

value *= 1 / voxel_volume  # eV m-3 per neutron
value *= source_strength  # eV m-3 s-1
value *= 1.602e-19  # J s-1 m-3 = W m-3

value = rmp.reshape_values_to_mesh_shape(heat_generation_mesh, value)
value[np.where(value < 0.01)] = np.nan

threshold = 3e8
value[np.where(value < threshold)] = np.nan

rmp.plot_regular_mesh_values(
    values=value,
    extent=rmp.get_tally_extent(heat_generation_mesh),
    rotate_plot=180,
    label="Heat generation (W m$^{-3}$)",
)

plt.savefig("heat_generation.png")

# depth distribution of helium generation
distribution = np.nanmean(value, axis=1)
distribution = distribution[np.logical_not(np.isnan(distribution))]
distribution *= 1e-6  # MW
height_mb = 2.5  # cm
x = np.linspace(0, height_mb, len(distribution), endpoint=True)
res = linregress(x[1:], distribution[::-1][1:])


plt.figure()
plt.scatter(x, distribution[::-1])
plt.plot(x, res.slope * x + res.intercept)
plt.xlabel("Distance from the top surface (cm)")
plt.ylabel("Heat generation  (MW m$^{-3}$)")
plt.ylim(bottom=0)
plt.savefig("heat_distribution.png")
