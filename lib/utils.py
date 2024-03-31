import numpy as np
import trimesh


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


def show_mesh(mesh):
    scene = trimesh.scene.Scene()
    scene.add_geometry(mesh)
    scene.show(smooth=False)


def show_manifold(manifold):
    show_mesh(manifold2trimesh(manifold))


def export_manifold(manifold, filename):
    manifold2trimesh(manifold).export(filename)
