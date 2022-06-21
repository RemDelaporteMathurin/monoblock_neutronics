import cadquery as cq
import paramak as pk

# construction


class Monoblock:
    def __init__(
        self,
        thickness,
        height,
        width,
        cucrzr_inner_radius,
        cucrzr_thickness,
        cu_thickness,
        w_thickness,
        gap,
        hollow=True,
        location=(0, 0, 0),
        normal=(0, 0, 1),
        xDir=None,
    ) -> None:
        """_summary_
        Args:
            thickness (float): thickness of the monoblock (mm)
            height (float): height of the monoblock in the Y direction (mm)
            width (float): width of the monoblock in the X direction (mm)
            cucrzr_inner_radius (float): inner radius of the CuCrZr pipe (mm)
            cucrzr_thickness (float): thickness of the CuCrZr pipe (mm)
            cu_thickness (float): thickness of the Cu interlayer (mm)
            w_thickness (float): thickness of W above the Cu
                (in the middle of the MB) (mm)
            gap (float): Poloidal gap between two monoblocks (mm)
        """
        self.thickness = thickness
        self.height = height
        self.width = width
        self.cucrzr_inner_radius = cucrzr_inner_radius
        self.cucrzr_thickness = cucrzr_thickness
        self.cu_thickness = cu_thickness
        self.w_thickness = w_thickness
        self.gap = gap
        self.location = location
        self.normal = normal
        self.xDir = xDir
        self.hollow = hollow
        self.plane = cq.Plane(self.location, normal=self.normal, xDir=self.xDir)
        self.make_solid()

    def make_solid(self):

        water = cq.Workplane(self.plane).cylinder(
            self.thickness * 2,
            self.cucrzr_inner_radius,
        )

        cucrzr = (
            cq.Workplane(self.plane)
            .cylinder(
                self.thickness + self.gap,
                self.cucrzr_inner_radius + self.cucrzr_thickness,
            )
            .cut(water)
        )

        copper = cq.Workplane(self.plane).cylinder(
            self.thickness,
            self.cucrzr_inner_radius + self.cucrzr_thickness + self.cu_thickness,
        )
        if self.hollow:
            copper = copper.cut(cucrzr).cut(water)

        translation_factor = (
            -self.height / 2
            + self.cucrzr_inner_radius
            + self.cucrzr_thickness
            + self.cu_thickness
            + self.w_thickness
        )

        tungsten = (
            cq.Workplane(self.plane)
            .move(-translation_factor, 0)
            .box(self.height, self.width, self.thickness)
            .cut(copper)
            .cut(cucrzr)
            .cut(water)
        )

        self.tungsten = tungsten
        self.copper = copper
        self.cucrzr = cucrzr
        self.water = water


if __name__ == "__main__":
    my_mb = Monoblock(
        thickness=1.2,
        height=2.5,
        width=2.3,
        cucrzr_inner_radius=0.6,
        cucrzr_thickness=0.15,
        w_thickness=0.5,
        cu_thickness=0.1,
        gap=0.0,
        location=(100, 0, 0),
        normal=(-1, 0, 0),
        hollow=True,
    )

    tungsten = pk.Shape(name="tungsten")
    tungsten.solid = my_mb.tungsten

    copper = pk.Shape(name="copper")
    copper.solid = my_mb.copper

    cucrzr = pk.Shape(name="cucrzr")
    cucrzr.solid = my_mb.cucrzr

    water = pk.Shape(name="water")
    water.solid = my_mb.water

    my_monoblock_paramak = pk.Reactor([tungsten, copper, cucrzr, water])

    my_monoblock_paramak.export_dagmc_h5m("monoblock.h5m")
    cq.exporters.export(my_mb.tungsten, "w.stl")
    cq.exporters.export(my_mb.copper, "copper.stl")
    cq.exporters.export(my_mb.cucrzr, "cucrzr.stl")
