from manifold3d import *

from lib.utils import *

cube_count = 200
random_points = np.random.uniform(low=-1, high=1, size=(cube_count, 3))
normalized_points = random_points / np.linalg.norm(random_points, axis=1, keepdims=True) * 5
random_axes = np.random.uniform(low=-1, high=1, size=(cube_count, 3))
random_angles = np.random.uniform(low=0, high=180, size=(cube_count, 3))

base_box = Manifold.cube([1, 1, 1], True)
sphere = Manifold.sphere(1.2, 100)

for point, angle in zip(normalized_points, random_angles):
    print(point, angle)
    base_box += Manifold.cube([1, 1, 1], True).rotate(angle).translate(point)

layer_height = 35
layer_count = 4
cube_count = 10

cubes = []
for i, radius in enumerate(np.linspace(90, 120, layer_count)):
    for theta in np.linspace(0, 2 * np.pi - (2 * np.pi / cube_count), cube_count):
        phase = i * 0.3
        x = radius * np.sin(theta + phase)
        y = radius * np.cos(theta + phase)
        cube = Manifold.cube([2, 2, 2], True)
        scale = np.random.uniform(15, 25, 1)[0]
        cube = cube.scale([scale] * 3)
        cube = cube.rotate(np.random.uniform(0, 90, 3))
        cubes.append(cube.translate([x, y, i * layer_height]))

cube_pot = sum(cubes[1:], cubes[0])

if True:
    set_min_circular_angle(5)

cube_pot += Manifold.cylinder(150, 80, 110).translate([0, 0, -40])
cube_pot -= Manifold.cylinder(180, 75, 110).translate([0, 0, -35])

export_manifold(cube_pot, "exported_stl/cube_topf.stl")
show_manifold(cube_pot)
