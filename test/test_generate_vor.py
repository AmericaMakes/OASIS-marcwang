import unittest

from OasisLib import GeoMesh, mesh_load, mesh_repair, mesh_save
from OasisLib import mesh_tetrahedralize, remesh_smooth
from OasisLib import MeshAttributesFlags, MeshIOFlags, polyhedral_mesher

class GeogramTest(unittest.TestCase):

    def setUp(self):
        self.stl_path = './test/test_artifact/Parameter_quality_nut_1.stl'
        self.out = './test/test_output/bench.obj'

    def test_generate_voronoi(self):
        m_in = GeoMesh()
        load_status = mesh_load(self.stl_path, m_in)
        self.assertTrue(load_status)

        mesh_repair(m_in)

        m_out = GeoMesh()
        remesh_smooth(m_in, m_out)

        tetra_status = mesh_tetrahedralize(m_out)
        self.assertTrue(tetra_status)

        voronoi = GeoMesh()
        status_vor = polyhedral_mesher(m_out, voronoi, tessallate_non_convex = True)
        self.assertTrue(status_vor)
        
        out_attr = MeshIOFlags()
        out_attr.set_attributes(MeshAttributesFlags.MESH_ALL_ATTRIBUTES)
        mesh_save(voronoi, self.out, out_attr)

    def test_access_attr(self):
        m_in = GeoMesh()
        load_status = mesh_load(self.stl_path, m_in)
        array = m_in.vertices.point(0)
        self.assertTrue(load_status)
        self.assertTrue(len(array) > 0)

if __name__ == '__main__':
    unittest.main()
