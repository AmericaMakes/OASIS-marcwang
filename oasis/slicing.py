import networkx as nt
import OasisLib as osl
import numpy as np
import trimesh

from oasis.hatching import straight_hatch

from shapely.geometry import Polygon
from typing import Dict, List
from scipy.spatial import KDTree
import random

def slice_object(filename: str, layer_thickness: float, h_width):

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

    contour_ls = []
    hatch_ls = []
    prev_poly = {}
    prev_kdt = None
    prev_g = None
    mhsl = osl.MeshHeightSlicer(voronoi)
    for i in np.arange(layer_thickness, bbs[2] - layer_thickness, layer_thickness):
        single = osl.GeoMesh()
        mhsl.get_layer(single, i)
        nb_f = single.facets.nb()
        current_poly = [convert_mesh_poly(single, f) for f in range(nb_f)]

        g, kd = init_g(current_poly, prev_poly, prev_kdt, prev_g, single)

        hatch_ls.append(process_hatch(g, current_poly, h_width))
        trimesh_2d_sl = trimesh_inst.section(plane_origin=[0, 0, i],
                                             plane_normal=[0, 0, 1])

        slice_2D, _ = trimesh_2d_sl.to_planar()
        contour_ls.append(slice_2D)

        prev_poly = current_poly
        prev_kdt = kd
        prev_g = g
    
    return hatch_ls, contour_ls

def init_g(current_poly: List[Polygon], prev_poly: List[Polygon],
           prev_kdt: KDTree, prev_g: nt.Graph, mesh: osl.GeoMesh):
    g = nt.Graph()
    data = []
    for f, poly in enumerate(current_poly):
        centroid = np.array(poly.centroid.coords).squeeze()
        data.append(centroid)
        w, ang = compute_prev_cost_and_angle(centroid, prev_kdt, prev_poly, prev_g)
        g.add_node(f, weight=w, angle = ang)

    t = KDTree(data)

    nb_f = mesh.facets.nb()
    for f in range(nb_f):
        c_range = range(mesh.facets.corners_begin(f),
                        mesh.facets.corners_end(f))
        for c in c_range:
            f_adj = mesh.facet_corners.adjacent_facet(c)

            if(f_adj != -1 and f_adj != 4294967295):
                g.add_edge(f, f_adj)
    return g, t


def compute_prev_cost_and_angle(c_centroid: np.ndarray, prev_kdt: KDTree,
                      prev_poly: List[Polygon], prev_g: nt.Graph):
    if prev_kdt is None or prev_g is None:
        return 0, random.uniform(0, 360)

    _, f_cl = prev_kdt.query(c_centroid)
    cost = prev_g.nodes[f_cl]['weight']
    angle = prev_g.nodes[f_cl]['angle'] + 67

    for f_adj in prev_g.adj[f_cl]:
        cost += prev_g.nodes[f_adj]['weight'] * 0.5
    return cost, angle


def convert_mesh_poly(mesh: osl.GeoMesh, f):
    c_it = range(mesh.facets.corners_begin(f), mesh.facets.corners_end(f))

    ring = []
    for c in c_it:
        p_1_id = mesh.facet_corners.vertex(c)
        p1 = mesh.vertices.point(p_1_id)
        ring.append(p1[0:2])

    return Polygon(ring)


def process_hatch(g: nt.Graph, current_poly: List[Polygon], hatch : float):
    # TODO sort on subset and pop to avoid reprocessing same
    n_path = len(current_poly)
    mod_g = g.copy()
    sorted_n = sorted(list(mod_g.nodes(data='weight')), key = lambda x : x[1])
    min_g = sorted_n[0]
    count = 0
    hatch_order = []
    while(count < n_path):
        poly = current_poly[min_g[0]]
        angle = mod_g.nodes[min_g]['angle']
        hatch = straight_hatch(poly, hatch, angle)
        hatch_order.append(hatch)

        c_area = poly.area

        for c in g.adj[min_g[0]].keys():
            g.nodes[c]['weight'] += c_area * 0.5

        if len(mod_g) > 1:
            mod_g.remove_node(min_g[0])
            sorted_n = sorted(list(mod_g.nodes(data='weight')), key = lambda x : x[1])
            min_g = sorted_n[0]
        
        count += 1

    return hatch_order
