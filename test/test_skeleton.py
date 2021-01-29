import unittest
import trimesh

from oasis.skeleton import Skeletonization, edge_collapse, vertex_collision, compute_edge_velocity

import numpy as np
import time


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
        skel.show_velocity_inner(0,0.1)

    def test_edge_collapse(self):
        tri = np.array([
            [0, 0],
            [2, 5],
            [4, 0]
        ], dtype=float)

        vel = np.array([
            [1, 0],
            [0, 0],
            [-1, 0]], dtype=float
        )

        test = edge_collapse(tri, vel)
        self.assertEqual(test[0], 2.0)
        self.assertEqual(test[1][0], 2)
        self.assertEqual(test[1][1], 0)
    
    def test_vertex_split_collision(self):
        tri = np.array([
            [0, 0],
            [2, 5],
            [4, 0]
        ], dtype=float)

        edges = np.roll(tri, shift = -1, axis = 0) - tri 
        wavefront = np.array([False, False, True])
        vel = np.array([
            [0, 1],
            [0,-1],
            [0, 1]], dtype=float
        )

        test = vertex_collision(tri, vel, edges, wavefront)

        self.assertEqual(test[0], 2.5)
        self.assertEqual(test[1], 1)
        self.assertEqual(test[2], 1)
    
    def test_vertex_flip_collision(self):
        tri = np.array([
            [0, 0],
            [2, 5],
            [4, 0]
        ], dtype=float)

        edges = np.roll(tri, shift = -1, axis = 0) - tri 
        wavefront = np.array([False, False, False])

        vel = np.array([
            [0, 1],
            [0,-1],
            [0, 1]], dtype=float
        )

        test = vertex_collision(tri, vel, edges, wavefront)

        self.assertEqual(test[0], 5.0)
        self.assertEqual(test[1], 1)
        self.assertEqual(test[2], 2)

    def test_compute_wavefront(self):
        skel = Skeletonization(self.polygons_2d[0])
        skel.compute_wavefront()


if __name__ == '__main__':
    unittest.main()
