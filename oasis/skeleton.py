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


#@jit(nopython=True, parallel=True)
def vertex_collision(tri: np.array, vel: np.array):
    edges = np.vstack(
        (tri[1] - tri[0],
         tri[2] - tri[1],
         tri[0] - tri[2])
    )

    signed_area = 0.5*LA.det(np.concatenate((tri.T, ((0, 0, 0))), axis = 0))
    edge_split = []
    for i in prange(3):
        if i < 2:
            vert_idx = i + 2
        else:
            vert_idx = 2 - i

        vert_vel = vel[vert_idx]
        if signed_area > 0:
            z_axis = np.array((0, 0, 1))
        else:
            z_axis = np.array((0, 0, -1))

        edge_motion = np.cross(z_axis,
            np.concatenate((edges[i],[0])))[:3]
        normed_edge_motion = edge_motion/LA.norm(edge_motion)
        normed_edge = edges[i]/LA.norm(edges[i])

        change_basis = LA.inv(np.hstack((normed_edge_motion, normed_edge)))
        vel_relative = np.dot(change_basis, vel[vert_idx])
        if vel_relative[0] < 0:
            distance = 2*np.abs(signed_area)/LA.norm(edges[i])
            time = np.abs(distance/vel_relative[0])
            edge_split.append((time, vert_idx))

    return edge_split

def compute_vertex_velocity(vert : np.array, holes : bool = False):
    edges = vert[1:,...] - vert[0:-1,...]

    if holes:
        dir = np.array((0,0,1))
    else:
        dir = np.array((0,0,-1))
    
    perp_dir = np.cross(edges, dir)[...,0:2]
    perp_dir = perp_dir/LA.norm(perp_dir, axis=-1)[...,np.newaxis]
    
    return perp_dir[0:,...] + np.roll(perp_dir, shift = 1, axis = 0)

class Skeletonization():
    def __init__(self, polygon: Polygon) -> None:
        precision = 1e-3
        # remove colinear segment
        outer = merge_colinear(polygon.exterior.coords, scale=precision)
        inners = [merge_colinear(hole.coords, scale=precision)
                  for hole in polygon.interiors]
        simplified_poly = Polygon(outer, holes=inners)

        self.poly = orient(simplified_poly)

        self.shell_bisector = compute_vertex_velocity(np.array(self.poly.exterior.coords))
        self.holes_bisector = [compute_vertex_velocity(
            np.array(hole.coords), holes=True) for hole in self.poly.interiors]

    def show_velocity_outer(self, scale = 1.0):
        polygon = load_path(self.poly)
        vert_1 = np.array(self.poly.exterior.coords)[:-1]
        vert_2 = vert_1 + self.shell_bisector*scale
        lines = np.nan_to_num(np.stack((vert_1, vert_2), axis=1))
        shapely_lines = load_path(lines)

        (shapely_lines + polygon).show()

    def show_velocity_inner(self, idx: int = 0, scale = 1.0):
        polygon = load_path(self.poly)
        vert_1 = np.array(self.poly.interiors[idx].coords)[:-1]
        vert_2 = vert_1 + self.holes_bisector[idx]*scale
        lines = np.nan_to_num(np.stack((vert_1, vert_2), axis=1))
        shapely_lines = load_path(lines)
        (shapely_lines + polygon).show()

    def compute_wavefront(self, step: float = 0.1):
        self.mesh_v, self.mesh_f = triangulate_polygon(
            self.poly, engine='earcut')

        velocity = np.concatenate(
            [self.shell_bisector, [self.shell_bisector[0]]])

        for hole in self.holes_bisector:
            velocity = np.concatenate([velocity, hole, [hole[0]]])

        for tri in self.mesh_f:
            vert = self.mesh_v[tri]
            vel = velocity[tri]
            print('test')
            #collapse_time = self.compute_collapse(vert, vel)

    def show_mesh(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(self.mesh_v[:, 0], self.mesh_v[:, 1], np.zeros(
            self.mesh_v.shape[0]), triangles=self.mesh_f)
        plt.show()
