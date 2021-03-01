import networkx as nt
import OasisLib as osl
import numpy as np
import trimesh

from oasis.hatching import straight_hatch

from shapely.geometry import Polygon
from typing import Dict

def slice_object(filename: str, layer_thickness: float):

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
    contour_ls = []
    hatch_ls = []
    prev_cell = {}
    prev_slice = None
    mhsl = osl.MeshHeightSlicer(voronoi)
    for i in np.arange(0, bbs[2], layer_thickness):
        single = osl.GeoMesh()
        f2c = mhsl.get_layer(single, i)
        c2f = dict((v,k) for k,v in f2c.items())

        g = update_graph(g, f2c, single)
        g = initialize_weight(g, prev_slice, f2c, prev_cell)

        hatch_ls.append(process_hatch(g, f2c, c2f, single))

        trimesh_2d_sl = trimesh_inst.section(plane_origin=[0, 0, i],
                                             plane_normal=[0, 0, 1])

        slice_2D, _ = trimesh_2d_sl.to_planar()
        contour_ls.append(slice_2D)
        
        prev_cell = f2c
        prev_slice = single

def convert_mesh_poly(mesh : osl.GeoMesh, f):
    c_it = range(mesh.facets.corner_begin(f), mesh.facets.corner_end(f))

    ring = []
    for c in c_it:
        p_1_id = mesh.facet_corners.vertex(c)
        p1 = mesh.vertices.point(p_1_id)
        ring.append(p1[0:2])

    return Polygon(ring)

def process_hatch(g : nt.Graph, f2c: Dict[int, int], c2f: Dict[int, int],  mesh: osl.GeoMesh):
    n_path = len(f2c)
    # TODO sort on subset and pop to avoid reprocessing same
    sorted_n = sorted(g.nodes(data=True), key = lambda x : x[1]['weight'])
    min_g = sorted_n[0]
    count = 0
    hatch_order = []
    while( count < n_path):
        f_id = c2f[min_g[0]]
        poly = convert_mesh_poly(mesh, f_id)
        # TODO setting input later
        hatch = straight_hatch(poly, 0.1, 67)
        hatch_order.append(hatch)

        c_area = poly.area

        for c, _ in g.adj[min_g[0]]:
            g[c]['weight'] += c_area*0.5

        sorted_n = sorted(g.nodes(data=True), key = lambda x : x[1]['weight'])
        min_g = sorted_n[0]
        count += 1
    
    return hatch_order

def update_graph(g: nt.Graph, id_map: Dict[int, int], mesh: osl.GeoMesh):
    facets = mesh.facets
    for f, c in id_map.items():
        if not g.has_node(c):
            g.add_node(c, weight = 0)

            for i in range(facets.nb_vertices(f)):
                f_adj = facets.adjacent(f, i)
                if f_adj != -1:
                    g.add_edge(id_map[f_adj], c)
    return g

def compute_area_poly(mesh : osl.GeoMesh, f : int):
    c_it = range(mesh.facets.corner_begin(f), mesh.facets.corner_end(f))
    area_sum = 0

    for c in c_it:
        p_1_id = mesh.facet_corners.vertex(c)
        p1 = mesh.vertices.point(p_1_id)

        c_2 = mesh.facets.next_corner_around_facet(f, c)
        p_2_id = mesh.facet_corners.vertex(c_2)
        p2 = mesh.vertices.point(p_2_id)

        area_sum += p1[0]*p2[1] - p1[1]*p2[0]

    return area_sum/2

def initialize_weight(g : nt.Graph, prev_layer : osl.GeoMesh, current_cell: Dict[int, int], prev_cell : Dict[int, int]):
    if len(prev_cell) == 0 or prev_layer is None:
        return 
    
    areas = {}
    for p_f, p_c in prev_cell.items():
        area_f = compute_area_poly(prev_layer, p_f)
        areas[p_c] = area_f

    n_g = g.copy()
    for c in current_cell.values():
        w = g[c]['weight']
        for k, _ in g.adj[c].items():
            if k in prev_cell.keys():
                w += areas[k]*0.5

        n_g[c]['weight'] = w

    return n_g
