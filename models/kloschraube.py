from manifold3d import *
from lib.utils import *

e = 5 * 2 / 3 ** .5
print(e)
c1 = Manifold.cube((10, e, 5), True)
c2 = Manifold.cube((10, e, 5), True).rotate((0, 0, 60))
c3 = Manifold.cube((10, e, 5), True).rotate((0, 0, 120))

hex = c1 + c2 + c3
hex2 = Manifold.compose([c1, c2, c3]).hull().translate((20, 0, 0))

hex3 = Manifold.batch_boolean([c1, c2, c3], OpType.Add).translate((40, 0, 0))

z1 = Manifold.cylinder(5, 7, 7, 100).translate((0, 0, 0))
z2 = Manifold.cylinder(5, 7, 7, 100).translate((20, 0, 0))
z3 = Manifold.cylinder(5, 7, 7, 100).translate((40, 0, 0))

show_manifold(
    Manifold.compose([z1 - hex, z2 - hex2, z3 - hex3])
)
