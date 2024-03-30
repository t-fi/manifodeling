from manifold3d import Manifold
import trimesh
import numpy as np
from time import time


# Helper to convert a Manifold into a Trimesh
def manifold2trimesh(manifold):
    mesh = manifold.to_mesh()

    if mesh.vert_properties.shape[1] > 3:
        vertices = mesh.vert_properties[:, :3]
        colors = (mesh.vert_properties[:, 3:] * 255).astype(np.uint8)
    else:
        vertices = mesh.vert_properties
        colors = None

    return trimesh.Trimesh(
        vertices=vertices, faces=mesh.tri_verts, vertex_colors=colors
    )


# Helper to display interactive mesh preview with trimesh
def showMesh(mesh):
    t0 = time()
    scene = trimesh.Scene()

    print(f"Scene took {(time() - t0) * 1000:.1f}ms")
    t0 = time()
    scene.add_geometry(mesh)
    print(f"add_geometry took {(time() - t0) * 1000:.1f}ms")
    t0 = time()
    scene.add_geometry(trimesh.creation.axis())
    print(f"axis took {(time() - t0) * 1000:.1f}ms")
    t0 = time()
    scene.show()
    print(f"show took {(time() - t0) * 1000:.1f}ms")
    t0 = time()


def showMesh2(mesh):
    scene = trimesh.scene.Scene()
    scene.add_geometry(mesh)
    print("showing")
    scene.show(smooth=False)


def posColors(pos, _):
    return [-p + 0.5 for p in pos] + [1.0]


t0 = time()

cube_count = 20
random_points = np.random.uniform(low=-1, high=1, size=(cube_count, 3))
normalized_points = random_points / np.linalg.norm(random_points, axis=1, keepdims=True) * 1
random_axes = np.random.uniform(low=-1, high=1, size=(cube_count, 3))
random_angles = np.random.uniform(low=0, high=180, size=(cube_count, 3))

base_box = Manifold.cube([1, 1, 1], True)

for point, angle in zip(normalized_points, random_angles):
    print(point, angle)
    base_box += Manifold.cube([1, 1, 1], True).rotate(angle).translate(point)

# Main

print(f"Generation took {(time() - t0) * 1000:.1f}ms")

t0 = time()
mesh = manifold2trimesh(base_box)

print(f"Mesh took {(time() - t0) * 1000:.1f}ms")
showMesh2(mesh)
