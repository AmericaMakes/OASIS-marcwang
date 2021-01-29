import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA

from shapely.geometry.polygon import Polygon, orient
from shapely.coords import CoordinateSequence
from trimesh.creation import triangulate_polygon
from trimesh.path.exchange.load import load_path

from trimesh.path.simplify import merge_colinear
from scipy.optimize import minimize

from typing import NewType, Tuple, List

from numba import jit, prange


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


# @jit(nopython=True, parallel=True)
def vertex_collision(vert: np.array, vert_vel: np.array, edges: np.array, wavefront: np.array) -> Tuple[float, int, int]:
    signed_area = 0.5 * \
        LA.det(np.concatenate(
            (vert, np.array([1, 1, 1])[..., np.newaxis]), axis=1))

    edge_t = np.array([np.inf, np.inf, np.inf])
    event_type = np.array([0,0,0])
    for i in prange(3):
        vert_v_i = vert_vel[i, ...]
        adj_e = edges[i-1]
        if i < 2:
            j = i+1
            op_e = edges[j, ...]
        else:
            j = 2-i
            op_e = edges[j, ...]

        if wavefront[j]:
            if signed_area < 0:
                dir = np.array((0, 0, 1))
            else:
                dir = np.array((0, 0, -1))
            perp_dir = np.cross(op_e, dir)
            opp_e_vel = (perp_dir/LA.norm(perp_dir))[0:2]
        else:
            opp_e_vel = np.array((0, 0))

        dist = adj_e - np.dot(adj_e, op_e) / \
            np.dot(op_e, op_e)*op_e

        rel_speed = opp_e_vel - vert_v_i
        t = np.dot(dist, dist)/np.dot(rel_speed, dist)

        if t > 0 and wavefront[j]:
            edge_t[i] = t
            event_type[i] = 1 # split
        elif t > 0:
            edge_t[i] = t
            event_type[i] = 2 # flip

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


@jit(nopython=True)
def compute_vertex_velocity(edge_vel: np.array) -> np.array:
    # assume last is identical to first
    edge_vel = edge_vel[0:-1, ...] + edge_vel[1:, ...]
    # add first to last
    return np.concatenate((edge_vel, np.expand_dims(edge_vel[0], axis=0)), axis=0)


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
