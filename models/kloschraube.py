from manifold3d import *

from lib.utils import *

main_body = Manifold.cylinder(10, 7, 15, 100)

handle = Manifold.cylinder(10, 3, 3, 100).translate((12, 0, 0))
handle += Manifold.cylinder(10, 3, 3, 100).translate((-12, 0, 0))

handle = handle.hull()

main_body += handle

# construct hex hole
e = 5 * 2 / 3 ** .5
c1 = Manifold.cube((10, e, 8), True)
c2 = Manifold.cube((10, e, 8), True).rotate((0, 0, 60))
hex_hole = Manifold.compose([c1, c2]).hull()

screw_hole = Manifold.cylinder(10, 3.5, 3.5, 100)

kloschraube = main_body - screw_hole - hex_hole

export_manifold(kloschraube, "../exported_stl/kloschraube.stl")

show_manifold(kloschraube)
