from manifold3d import *

from lib.utils import *

main_body = Manifold.cylinder(176, 25, 25, 100)
main_body += Manifold.cylinder(170, 19.9, 19.9, 100).translate((0, 0, 10))

main_body += Manifold.cylinder(6, 89, 89, 100, center=True).rotate((90, 0, 0)).translate((0, 0, 89))
main_body -= Manifold.cylinder(60, 83, 83, 100, center=True).rotate((90, 0, 0)).translate((0, 0, 89))

main_body -= Manifold.cylinder(50, 20.1, 20.1, 100)
main_body -= Manifold.cylinder(500, 15, 15, 100)

export_manifold(main_body, "../exported_stl/lampenspacer.stl")
show_manifold(main_body)
