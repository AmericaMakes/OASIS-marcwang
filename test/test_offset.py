import unittest
import trimesh

import oasis.clipper_offsetting as clo

import matplotlib.pyplot as plt

from trimesh.path import polygons

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
    

    def test_offsetting(self):

        clo_obj = clo.ClipperOffsetting(self.polygons_2d[0])

        offsets = clo_obj.get_offset( 1.0, 3)
        outer_path = trimesh.load_path(self.polygons_2d[0])
        #outer_path.show()

        ls_path = [trimesh.load_path(poly) for poly in offsets]

        for p in ls_path:
            outer_path += p

        outer_path.show()

        