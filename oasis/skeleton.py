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


@jit(nopython=True)
def edge_collapse(tri: np.array, vel: np.array) -> List[Tuple[float, np.array]]:
    edge_collapse = []
    pair_edges = [(0, 1), (1, 2), (2, 0)]
    for i in range(3):
        p = pair_edges[i]
        min_t = min_edge_length_time(
            tri[p[0]], tri[p[1]], vel[p[0]], vel[p[1]])
        if min_t > 0:
            length = compute_edge_length(
                min_t, tri[p[0]], tri[p[1]], vel[p[0]], vel[p[1]])
            if length < 1e-3:
                edge_collapse.append((min_t, p))
    return edge_collapse


# @jit(nopython=True, parallel=True)
def vertex_collision(tri_vert: np.array, vert_vel: np.array, tri_edges: np.array, edges_vel: np.array):

    signed_area = 0.5 * \
        LA.det(np.concatenate(
            (tri_vert, np.array([0, 0, 0])[..., np.newaxis]), axis=1))

    edge_split = []
    for i in range(3):
        vert_i = tri_vert[i, ...]
        if i < 2:
            opposite_edge = tri_edges[i+1, ...]
            opposite_edge_vel = edges_vel[i+1, ...]
        else:
            opposite_edge = tri_edges[2-i, ...]
            opposite_edge_vel = edges_vel[2-1, ...]

    return edge_split


def compute_edge_velocity(vert: np.array, holes: bool = False):
    edges = vert[1:, ...] - vert[0:-1, ...]

    if holes:
        dir = np.array((0, 0, 1))
    else:
        dir = np.array((0, 0, -1))

    perp_dir = np.cross(edges, dir)[..., 0:2]
    perp_dir = perp_dir/LA.norm(perp_dir, axis=-1)[..., np.newaxis]
    return perp_dir


def compute_vertex_velocity(edge_vel: np.array):
    return edge_vel[0:, ...] + np.roll(edge_vel, shift=1, axis=0)


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
        vert_2 = vert_1 + self.shell_vert_vel*scale
        lines = np.nan_to_num(np.stack((vert_1, vert_2), axis=1))
        shapely_lines = load_path(lines)

        (shapely_lines + polygon).show()

    def show_velocity_inner(self, idx: int = 0, scale=1.0):
        polygon = load_path(self.poly)
        vert_1 = np.array(self.poly.interiors[idx].coords)[:-1]
        vert_2 = vert_1 + self.holes_vert_vel[idx]*scale
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
