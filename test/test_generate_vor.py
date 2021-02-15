import unittest

from OasisLib import Mesh, mesh_load, mesh_repair

class SkeletonizationTest(unittest.TestCase):

    def setUp(self):
        stl_path = './test/test_artifact/3DBenchy.stl'
        self.m = Mesh()
        mesh_load(stl_path, self.m)

    def test_pt_sampler(self):
        mesh_repair(self.m)
        
