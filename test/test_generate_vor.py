import unittest

from OasisLib import GeoMesh, mesh_load, mesh_repair, mesh_save
from OasisLib import mesh_tetrahedralize, remesh_smooth
from OasisLib import MeshAttributesFlags, MeshIOFlags, polyhedral_mesher
from OasisLib import MeshHeightSlicer

class GeogramTest(unittest.TestCase):

    def setUp(self):
        self.stl_path = './test/test_artifact/3DBenchy.stl'

        m_in = GeoMesh()
        load_status = mesh_load(self.stl_path, m_in)
        self.assertTrue(load_status)

        mesh_repair(m_in)

        self.tet = GeoMesh()
        remesh_smooth(m_in, self.tet, 10000)

        tetra_status = mesh_tetrahedralize(self.tet)
        self.assertTrue(tetra_status)

        self.voronoi = GeoMesh()
        status_vor = polyhedral_mesher(self.tet, self.voronoi, 1000,
            "tets_voronoi_boundary", 0.001, 5, 30, True, 0.0)
        self.assertTrue(status_vor)

        self.out_attr = MeshIOFlags()
        self.out_attr.set_attributes(MeshAttributesFlags.MESH_ALL_ATTRIBUTES)

    def test_generate_voronoi(self):
        out = './test/test_output/bench.obj'
        mesh_save(self.voronoi, out, self.out_attr)

    def test_single_layer(self):
        single = GeoMesh()
        self.assertTrue(self.voronoi.vertices.nb() > 0)
        mhsl = MeshHeightSlicer(self.voronoi)
        id_mapping = mhsl.get_layer(single, 20.0)

        out = './test/test_output/single_layer.obj'
        mesh_save(single, out, self.out_attr)

        self.assertTrue(len(id_mapping) > 0)
        self.assertTrue(single.vertices.nb() > 0)

if __name__ == '__main__':
    unittest.main()
