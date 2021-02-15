import unittest
import trimesh

from shapely.affinity import translate

import oasis.hatching as ohat
from oasis.util import clean_polygon

class HatchingTest(unittest.TestCase):

    def setUp(self):
        stl_path = './test/test_artifact/3DBenchy.stl'

        self.obj = trimesh.load(stl_path)

        height = (0, 0, 10)
        sec = self.obj.section(plane_origin=height,
                          plane_normal=[0, 0, 1])

        slice_2d, _ = sec.to_planar()

        slice_2d.merge_vertices()
        self.polygons_2d = [translate(poly, xoff=25, yoff=20) for poly in slice_2d.polygons_full]
    
    def test_pt_sampler(self):
        sampler = ohat.PointSampler(self.obj, 10)
        pt = sampler.generate_point(1000)
        pt_c = trimesh.points.PointCloud(pt, colors='blue')
        pt_c.show()
    
    def test_voronoi(self):
        v_infill = ohat.VoronoiInfill(self.obj, 10)
        v_infill.get_slice(10)
    
    def test_hatching(self):
        cleaned_poly = clean_polygon(self.polygons_2d[0])
        trimmed_hatch = ohat.straight_hatch(cleaned_poly, 0.1, 45)

        poly_path = trimesh.load_path(cleaned_poly)
        poly_path += trimesh.load_path(trimmed_hatch)
        poly_path.show()

if __name__ == '__main__':
    unittest.main()
