import unittest
import trimesh

from shapely.affinity import translate

import oasis.hatching as ohat
from oasis.util import clean_polygon

class SkeletonizationTest(unittest.TestCase):

    def setUp(self):
        stl_path = './test/test_artifact/3DBenchy.stl'

        obj = trimesh.load(stl_path)

        height = (0, 0, 10)
        sec = obj.section(plane_origin=height,
                          plane_normal=[0, 0, 1])

        slice_2d, _ = sec.to_planar()

        slice_2d.merge_vertices()
        self.polygons_2d = [translate(poly, xoff=25, yoff=20) for poly in slice_2d.polygons_full]
        
    
    def test_hatching(self):
        hat = ohat.Hatching(0)
        cleaned_poly = clean_polygon(self.polygons_2d[0])
        trimmed_hatch = hat.straight_hatch(cleaned_poly, 0.1, 45)

        poly_path = trimesh.load_path(cleaned_poly)
        poly_path += trimesh.load_path(trimmed_hatch)
        poly_path.show()