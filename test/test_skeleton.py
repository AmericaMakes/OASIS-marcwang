import unittest
import trimesh

from oasis.skeleton import Skeletonization

import numpy as np


class SkeletonizationTest(unittest.TestCase):

    def setUp(self):
        stl_path = './test/test_artifact/3DBenchy.stl'

        obj = trimesh.load(stl_path)

        height = (0, 0, 10)
        sec = obj.section(plane_origin=height,
                          plane_normal=[0, 0, 1])

        slice_2d, _ = sec.to_planar()

        self.polygons_2d = slice_2d.polygons_full

    def test_generate_bisector(self):
        skel = Skeletonization(self.polygons_2d[0])
        skel.show_velocity_outer()
        skel.show_velocity_inner(0)


if __name__ == '__main__':
    unittest.main()
