from manifold3d import Manifold
import numpy as np
from lib.utils import show_manifold

outer_wall = Manifold.cylinder(100, 100, 120)
inner_wall = Manifold.cylinder(95, 90, 110)

outer_wall -= inner_wall.translate((0,0,10))

show_manifold(outer_wall)
