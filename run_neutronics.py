import openmc
import openmc_dagmc_wrapper as odw
import openmc_plasma_source as ops

import neutronics_material_maker as nmm

import math

mb_height = 2.5
mb_thickness = 1.2
mb_width = 2.3
d = 100  # distance of the source from the monoblock top in cm

geometry = odw.Geometry(
    "monoblock.h5m",
    graveyard_box=(
        (-mb_thickness / 2, -mb_width / 2, -mb_height / 2),
        (mb_thickness / 2, mb_width / 2, mb_height / 2 + d),
    ),
    reflective_angles=(math.radians(-1), math.radians(1)),
)

water = nmm.Material(
    "H2O",
    density=0.9,
    chemical_equation="H2O",
    density_unit="g/cm3",
    percent_type="ao",
)

materials = odw.Materials(
    h5m_filename=geometry.h5m_filename,
    correspondence_dict={
        "copper": "copper",
        "cucrzr": "copper",  # TODO change this?
        "tungsten": "tungsten",
        "water": water,
    },
)


settings = odw.FusionSettings()
settings.batches = 50
settings.particles = 1000000

my_source = ops.FusionRingSource(
    fuel="DT",
    radius=100,
    z_placement=mb_height / 2 + d,
    angles=(math.radians(-1), math.radians(1)),
)
settings.source = my_source


tally = odw.CellTally(
    "(n,Xa)",
    target="tungsten",
    materials=materials,
    name="total_He_generation_in_tungsten",
)
tally2 = odw.CellTally(
    "(n,Xa)", target="copper", materials=materials, name="total_He_generation_in_copper"
)

regular_mesh_tally_helium = odw.MeshTally2D(
    "(n,Xa)",
    plane="yz",
    bounding_box=[
        (100 - mb_thickness / 2 * 1.1, -mb_width / 2 * 1.1, -mb_height / 2 * 1.1),
        (100 + mb_thickness / 2 * 1.1, mb_width / 2 * 1.1, mb_height / 2 * 1.1),
    ],
    plane_slice_location=[100 + mb_thickness / 2, 100 - mb_thickness / 2],
    mesh_resolution=(50, 50),
)
tungsten = [mat for mat in materials if mat.name == "tungsten"][0]
regular_mesh_tally_helium.filters.append(openmc.MaterialFilter(tungsten))


regular_mesh_tally_heating = odw.MeshTally2D(
    "heating",
    plane="yz",
    bounding_box=[
        (100 - mb_thickness / 2 * 1.1, -mb_width / 2 * 1.1, -mb_height / 2 * 1.1),
        (100 + mb_thickness / 2 * 1.1, mb_width / 2 * 1.1, mb_height / 2 * 1.1),
    ],
    plane_slice_location=[100 + mb_thickness / 2, 100 - mb_thickness / 2],
    mesh_resolution=(50, 50),
)
regular_mesh_tally_heating.filters.append(openmc.MaterialFilter([tungsten]))

tallies = openmc.Tallies(
    [tally, tally2, regular_mesh_tally_helium, regular_mesh_tally_heating]
)

my_model = openmc.Model(
    materials=materials, geometry=geometry, settings=settings, tallies=tallies
)

if __name__ == "__main__":
    statepoint_file = my_model.run(tracks=False)
