import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA

from shapely.geometry.polygon import Polygon, orient
from trimesh.creation import triangulate_polygon
from trimesh.path.exchange.load import load_path

from trimesh.path.simplify import merge_colinear

from typing import Tuple, List

from numba import jit, prange

from enum import IntEnum
import heapq as hq

import networkx as nx


class TriEvent(IntEnum):
    No = 0
    Edge = 1
    Split = 2
    Flip = 3
    Multi = 4
    Removed = 5


@jit(nopython=True)
def compute_edge_length(t: float, v_1: np.array, v_2: np.array, vel_1: np.array, vel_2: np.array) -> float:
    return LA.norm((v_1 - v_2) + (vel_1 - vel_2)*t)


@jit(nopython=True)
def min_edge_length_time(v_1: np.array, v_2: np.array, vel_1: np.array, vel_2: np.array) -> float:
    t = -1 * np.sum(2*(v_1 - v_2) *
                    (vel_1 - vel_2))/np.sum(2*(vel_1 - vel_2)**2)
    return t


@jit(nopython=True, parallel=True)
def edge_collapse(tri: np.array, vel: np.array) -> Tuple[float, Tuple[int, int]]:
    edge_collapse = []
    pair_edges = [(0, 1), (1, 2), (2, 0)]
    edge_t = np.array([np.inf, np.inf, np.inf])
    for i in prange(3):
        p = pair_edges[i]
        min_t = min_edge_length_time(
            tri[p[0]], tri[p[1]], vel[p[0]], vel[p[1]])
        if min_t > 0:
            length = compute_edge_length(
                min_t, tri[p[0]], tri[p[1]], vel[p[0]], vel[p[1]])
            if length < 1e-3:
                edge_t[i] = min_t
                edge_collapse.append((min_t, p))

    min_idx = np.argmin(edge_t)
    return (edge_t[min_idx], pair_edges[min_idx])


@jit(nopython=True)
def vertex_collision(vert_vel: np.array, edges: np.array, wavefront: np.array) -> Tuple[float, int, int]:
    edge_t = np.array([np.inf, np.inf, np.inf])
    event_type = np.array([0, 0, 0])
    for i in range(3):
        vert_v_i = vert_vel[i, ...]
        adj_e = edges[i-1]

        j = (i + 1) % 3
        op_e = edges[j, ...]

        if wavefront[j]:
            # assume tri is always ccw ordered
            dir = np.array((0.0, 0.0, - 1.0))
            # perp points inward
            perp_dir = np.cross(op_e, dir)
            opp_e_vel = (perp_dir/LA.norm(perp_dir))[0:2]
        else:
            opp_e_vel = np.array((0.0, 0.0))

        dist = adj_e - np.dot(adj_e, op_e) / \
            np.dot(op_e, op_e)*op_e

        rel_speed = opp_e_vel - vert_v_i

        perp_speed = np.dot(rel_speed, dist)

        if perp_speed != 0:
            t = np.dot(dist, dist)/perp_speed
        else:
            t = 0

        if t > 0 and wavefront[j]:
            edge_t[i] = t
            event_type[i] = TriEvent.Split  # split
        elif t > 0:
            edge_t[i] = t
            event_type[i] = TriEvent.Flip  # flip

    min_idx = np.argmin(edge_t)
    return (edge_t[min_idx], min_idx, event_type[min_idx])


@jit(nopython=True, parallel=True)
def compute_edge_velocity(vert: np.array, holes: bool = False) -> np.array:
    # assume last is identical to first
    edges = vert[1:, ...] - vert[0:-1, ...]

    if holes:
        dir = np.array((0, 0, 1))
    else:
        dir = np.array((0, 0, -1))

    perp_dir = np.cross(edges, dir)[..., 0:2]
    normed_perp = np.zeros_like(perp_dir)
    for i in prange(normed_perp.shape[0]):
        normed_perp[i, ...] = perp_dir[i, ...]/LA.norm(perp_dir[i, ...])
    # add first to last
    return np.concatenate((normed_perp, np.expand_dims(normed_perp[0], axis=0)), axis=0)


# @jit(nopython=True)
def compute_vertex_velocity(edge_vel: np.array) -> np.array:
    # assume last is identical to first
    edge_vel = np.roll(edge_vel[0:-1, ...] +
                       edge_vel[1:, ...], shift=1, axis=0)
    # add first to last
    return np.concatenate((edge_vel, np.expand_dims(edge_vel[0], axis=0)), axis=0)

def _check_if_node_last(g : nx.DiGraph, v1_idx : int, v2_idx : int, poly_loop: List[Tuple[int, int, int]]):
    if not g.has_node(v1_idx):
        l_id = g.nodes[v2_idx]['l_id']
        return poly_loop[l_id][0]
    else:
        return v1_idx

def create_graph(mesh_v: np.array, vel_v: np.array, mesh_f: np.array, poly_loop: List[Tuple[int, int, int]]):
    g = nx.DiGraph()
    m_g = nx.Graph()
    for l_id, loop in enumerate(poly_loop):
        interval = range(loop[0], loop[1])
        hole = loop[2]
        ls_nodes = [(i, {'coord': mesh_v[i, ...],
                         'velocity': vel_v[i, ...], 'l_id': l_id}) for i in interval]
        g.add_nodes_from(ls_nodes)

        ls_edges = [(i, (i+1) % loop[1] + loop[0]*(int((i+1)/loop[1])),
                     {'hole': hole, 'wave': True, 'l_id': l_id}) for i in interval]
        g.add_edges_from(ls_edges)

    for t_idx, tri in enumerate(mesh_f):
        m_g.add_node(t_idx)
        n_tri_v = np.zeros(3, dtype = int)
        for i in range(3):
            v_1_id = tri[i]
            v_2_id = tri[(i+1) % 3]

            v_1_id = _check_if_node_last(g, v_1_id, v_2_id, poly_loop)
            v_2_id = _check_if_node_last(g, v_2_id, v_1_id, poly_loop)

            if g.has_edge(v_1_id, v_2_id):
                g[v_1_id][v_2_id]['t_id'] = t_idx
            else:
                g.add_edge(v_1_id, v_2_id)
                g[v_1_id][v_2_id]['wave'] = False
                g[v_1_id][v_2_id]['t_id'] = t_idx

            if g.has_edge(v_2_id, v_1_id):
                if not g[v_2_id][v_1_id]['wave']:
                    adj_t_idx = g[v_2_id][v_1_id]['t_id']
                    m_g.add_edge(t_idx, adj_t_idx)
            
            n_tri_v[i] = v_1_id
        m_g.nodes[t_idx]['vertex'] = n_tri_v

    return g, m_g

def compute_wavefront(g_v : nx.DiGraph, g_m : nx.Graph):
    collapse_heap = []
    event_list = []

    for tri_idx in g_m.nodes:
        v_idx = g_m[tri_idx]['vertex']
        tri_v = np.array([g_v[v]['coord'] for v in v_idx]) 
        tri_v_vel = np.array([g_v[v]['velocity'] for v in v_idx])

        tri_edge = np.roll(tri_v, shift=-1, axis=0) - tri_v
        wavefront = np.array([g_v[v_idx[i]][v_idx[(i + 1) % 3]]['wave'] for i in range(3)])

        e_collapse = edge_collapse(tri_v, tri_v_vel)
        v_collision = vertex_collision(tri_v_vel, tri_edge, wavefront)

        if e_collapse[0] < v_collision[0]:
            v_pair_id = e_collapse[1]
            hq.heappush(collapse_heap, (e_collapse[0], tri_idx))
            event_list.append((TriEvent.Edge, v_pair_id))
        elif e_collapse[0] > v_collision[0]:
            v_id = v_idx[v_collision[1]]
            hq.heappush(collapse_heap, (v_collision[0], tri_idx))
            event_list.append((v_collision[2], v_id))
        else:
            if e_collapse[0] == np.inf and v_collision[0] == np.inf:
                hq.heappush(collapse_heap, (np.inf, tri_idx))
                event_list.append((TriEvent.No,))
            else:
                edge_id = e_collapse[1]
                v_id = v_idx[v_collision[1]]
                hq.heappush(collapse_heap, (e_collapse[0], tri_idx))
                event_list.append(
                    (TriEvent.Multi, ((TriEvent.Edge, edge_id), (v_collision[2], v_id))))

    processed_event = []  # (tri left, t collapse, sorted tri_idx)

    while len(collapse_heap) > 0:
        elem = hq.heappop(collapse_heap)
        event = event_list[elem[1]]
        tri_idx = g_m[elem[1]]['vertex']

        if event[0] == TriEvent.Edge:
            v_pair = event[1]
            g_v = nx.contracted_edge(g_v, v_pair)
            
        elif event[0] == TriEvent.Split:
            pass

    return collapse_heap


class Skeletonization():
    def __init__(self, polygon: Polygon) -> None:
        precision = 1e-3
        # remove colinear segment
        outer = merge_colinear(polygon.exterior.coords, scale=precision)
        inners = [merge_colinear(hole.coords, scale=precision)
                  for hole in polygon.interiors]
        simplified_poly = Polygon(outer, holes=inners)

        self.poly = orient(simplified_poly)

        self.shell_edge_vel = compute_edge_velocity(
            np.array(self.poly.exterior.coords))
        self.holes_edge_vel = [compute_edge_velocity(
            np.array(hole.coords), holes=True) for hole in self.poly.interiors]

        self.shell_vert_vel = compute_vertex_velocity(self.shell_edge_vel)
        self.holes_vert_vel = [compute_vertex_velocity(
            h_edge_vel) for h_edge_vel in self.holes_edge_vel]

    def show_velocity_outer(self, scale=1.0):
        polygon = load_path(self.poly)
        vert_1 = np.array(self.poly.exterior.coords)[:-1]
        vert_2 = vert_1 + self.shell_vert_vel[:-1]*scale
        lines = np.nan_to_num(np.stack((vert_1, vert_2), axis=1))
        shapely_lines = load_path(lines)

        (shapely_lines + polygon).show()

    def show_velocity_inner(self, idx: int = 0, scale=1.0):
        polygon = load_path(self.poly)
        vert_1 = np.array(self.poly.interiors[idx].coords)[:-1]
        vert_2 = vert_1 + self.holes_vert_vel[idx][:-1]*scale
        lines = np.nan_to_num(np.stack((vert_1, vert_2), axis=1))
        shapely_lines = load_path(lines)
        (shapely_lines + polygon).show()

    def compute_wavefront(self, step: float = 0.1):
        self.mesh_v, self.mesh_f = triangulate_polygon(
            self.poly, engine='earcut')

    def show_mesh(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(self.mesh_v[:, 0], self.mesh_v[:, 1], np.zeros(
            self.mesh_v.shape[0]), triangles=self.mesh_f)
        plt.show()
