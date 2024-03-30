from sdf import *
import numpy as np

random_points = np.random.uniform(low=-10, high=10, size=(10,3))
random_axes = np.random.uniform(low=-1, high=1, size=(10,3))
random_angles = np.random.uniform(low=0, high=np.pi, size=(10))
normalized_axes = random_axes / np.linalg.norm(random_points, axis=1, keepdims=True)

base_box = box(1.5)

for point, angle, axis in zip(random_points, random_angles, normalized_axes):
    print(point, angle, axis)
    base_box |= box(1.5).rotate(angle, axis).translate(point)

base_box.save('out.stl', step=0.05)