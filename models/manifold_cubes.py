from manifold3d import Manifold
import numpy as np
from lib.utils import show_manifold

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

#base_box ^= sphere

# Main

show_manifold(base_box)