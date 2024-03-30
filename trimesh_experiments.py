import numpy as np
import trimesh
colors = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
]
vertices = colors
scene = trimesh.scene.Scene()
for v, c in zip(vertices, colors):
    sphere = trimesh.primitives.Sphere(center=v, radius=0.1)
    # add point clouds - primitives don't render well
    vertex_colors = np.tile(
        np.expand_dims(c, axis=0)*255, (sphere.vertices.shape[0], 1))
    sphere.visual.vertex_colors = vertex_colors
    scene.add_geometry(sphere)

scene.show()