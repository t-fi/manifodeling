from manifold3d import *
from lib.utils import show_manifold, export_manifold

if True:
    set_min_circular_angle(5)


def bowl(height, radius_low, radius_high, thickness):
    outer_wall = Manifold.cylinder(height, radius_low, radius_high)
    inner_wall = Manifold.cylinder(height - thickness, radius_low - thickness, radius_high - thickness)
    return outer_wall - inner_wall.translate((0, 0, thickness))


lower_bowl = bowl(20, 90, 100, 2)

export_manifold(lower_bowl, filename="../exported_stl/lower_bowl.stl")
show_manifold(lower_bowl)
