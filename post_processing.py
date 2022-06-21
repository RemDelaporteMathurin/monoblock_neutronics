import openmc
import regular_mesh_plotter as rmp
from openmc_mesh_tally_to_vtk import write_mesh_tally_to_vtk
import matplotlib.pyplot as plt
import openmc_tally_unit_converter as otuc

result = openmc.StatePoint("statepoint.20.h5")
tally1 = result.get_tally(name="tungsten_(n,Xa)")

print(tally1.mean)
print(tally1.std_dev)

tally2 = result.get_tally(name="(n,Xa)_on_2D_mesh_yz")

value = tally2.mean

voxel_volume = otuc.compute_volume_of_voxels(tally2)  # cm3
print(voxel_volume)

voxel_volume *= 1e-6  # m3
print(voxel_volume)

source_strength = otuc.find_source_strength(
    fusion_energy_per_second_or_per_pulse=1000e6
)  # n/s
print(source_strength)


value *= 1 / voxel_volume  # m-3 per neutron
value *= source_strength  # m-3 s-1

value = rmp.reshape_values_to_mesh_shape(tally2, value)

rmp.plot_regular_mesh_values(
    values=value,
    extent=rmp.get_tally_extent(tally2),
    rotate_plot=180,
)

plt.savefig("openmc_mesh_tally_plot.png")
