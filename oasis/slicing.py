import networkx as nt
import OasisLib as osl 
import numpy as np
import trimesh
from typing import Dict

def slice_object(filename : str, layer_thickness : float):

    trimesh_inst = trimesh.load(filename)

    m_in = osl.GeoMesh()
    load_status = osl.mesh_load(filename, m_in)
    assert(load_status)

    osl.mesh_repair(m_in)
    tet = osl.GeoMesh()
    osl.remesh_smooth(m_in, tet, 10000)
    tetra_status = osl.mesh_tetrahedralize(tet)
    assert(tetra_status)

    voronoi = osl.GeoMesh()
    status_vor = osl.polyhedral_mesher(tet, voronoi, 1000,
            "tets_voronoi_boundary", 0.001, 5, 30, True, 0.0)
    assert(status_vor)
    bbs = trimesh_inst.extents

    g = nt.Graph()
    for i in np.arange(0, bbs[2], layer_thickness):
        single = osl.GeoMesh()
        mhsl = osl.MeshHeightSlicer(voronoi)
        id_mapping = mhsl.get_layer(single, i)

        build_graph(g, id_mapping, single)

        trimesh_2d_sl = trimesh_inst.section(plane_origin=[0,0,i], 
                     plane_normal=[0,0,1])
        
        slice_2D, _ = trimesh_2d_sl.to_planar()



def build_graph( g : nt.Graph, id_map : Dict[int,int], mesh : osl.GeoMesh):
    pass

