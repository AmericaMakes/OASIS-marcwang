import unittest
import trimesh

from oasis.skeleton import Skeletonization, compute_edge_velocity, compute_vertex_velocity, edge_collapse, vertex_collision, compute_wavefront
from trimesh.creation import triangulate_polygon
from shapely.geometry import Polygon
import numpy as np

import matplotlib.pyplot as plt
import geopandas as gpd


def plot_mesh(mesh_v, mesh_f):
    fig, ax = plt.subplots(figsize=(10, 10))

    def check_duplicate(ls, item):
        in_list = False
        for j in ls:
            if np.all(j == item):
                in_list = True
        return in_list

    added_v = []
    for i, v in enumerate(mesh_v):
        if not check_duplicate(added_v, v):
            ax.annotate(i, xy=v, color='b')
            added_v.append(v)
        else:
            ax.annotate(
                "%d*" % i, xy=v, xytext=(v[0]-0.2, v[1]-0.2), color='r', arrowprops={'arrowstyle': 'simple'})

    for i, tri in enumerate(mesh_f):
        vert_v = Polygon(mesh_v[tri])
        p = gpd.GeoSeries(vert_v)
        p.plot(ax=ax)

        centroid = vert_v.representative_point()

        ax.annotate("%d \n (%d,%d,%d)" %
                    (i, tri[0], tri[1], tri[2]), xy=(centroid.x, centroid.y))

    return fig, ax


class SkeletonizationTest(unittest.TestCase):

    def setUp(self):
        stl_path = './test/test_artifact/3DBenchy.stl'

        obj = trimesh.load(stl_path)

        height = (0, 0, 10)
        sec = obj.section(plane_origin=height,
                          plane_normal=[0, 0, 1])

        slice_2d, _ = sec.to_planar()

        slice_2d.merge_vertices()
        self.polygons_2d = slice_2d.polygons_full

    def test_generate_bisector(self):
        skel = Skeletonization(self.polygons_2d[0])
        skel.show_velocity_outer(0.5)
        skel.show_velocity_inner(0, 0.1)

    def test_edge_collapse(self):
        tri = np.array([
            [0, 0],
            [4, 0],
            [2, 5]
        ], dtype=float)

        vel = np.array([
            [1, 0],
            [-1, 0],
            [0, 0]], dtype=float
        )

        test = edge_collapse(tri, vel)
        self.assertEqual(test[0], 2.0)
        self.assertEqual(test[1][0], 0)
        self.assertEqual(test[1][1], 1)

    def test_vertex_split_collision(self):
        tri = np.array([
            [0, 0],
            [4, 0],
            [2, 5]
        ], dtype=float)

        edges = np.roll(tri, shift=-1, axis=0) - tri
        wavefront = np.array([False, False, True])
        vel = np.array([
            [0, 1],
            [0, 1],
            [0, -1]], dtype=float
        )

        test = vertex_collision(tri, vel, edges, wavefront)

        self.assertEqual(test[0], 2.5)
        self.assertEqual(test[1], 1)
        self.assertEqual(test[2], 1)

    def test_vertex_flip_collision(self):
        tri = np.array([
            [0, 0],
            [4, 0],
            [2, 5],
        ], dtype=float)

        edges = np.roll(tri, shift=-1, axis=0) - tri
        wavefront = np.array([False, False, False])

        vel = np.array([
            [0, 1],
            [0, 1],
            [0, -1]], dtype=float
        )

        test = vertex_collision(tri, vel, edges, wavefront)

        self.assertEqual(test[0], 5.0)
        self.assertEqual(test[1], 1)
        self.assertEqual(test[2], 2)

    def test_compute_wavefront(self):
        poly = Polygon([(0, 0), (0, 1), (1.5, 3.5), (5.0, 6.0),
                        (7.0, 0.5), (4.0, 0.0), (3.5, 1.5), (3.0, 0.0)])
        # plt.plot(*poly.exterior.xy)
        # plt.show()
        mesh_v, mesh_f = triangulate_polygon(
            poly, engine='earcut')

        poly_loop = [(0,9,0)]
        edge_vel = compute_edge_velocity(mesh_v)
        vert_vel = compute_vertex_velocity(edge_vel)

        plot_mesh(mesh_v, mesh_f)
        plt.show()
        compute_wavefront(mesh_v, vert_vel, mesh_f, poly_loop)


if __name__ == '__main__':
    unittest.main()
