import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA

from shapely.geometry.polygon import Polygon, orient
from shapely.coords import CoordinateSequence
from trimesh.creation import triangulate_polygon
from trimesh.path.exchange.load import load_path

from trimesh.path.simplify import merge_colinear


class Skeletonization():
    def __init__(self, polygon: Polygon) -> None:
        precision = 1e-3
        #remove colinear segment
        outer = merge_colinear(polygon.exterior.coords, scale=precision)
        inners = [merge_colinear(hole.coords, scale=precision)
                  for hole in polygon.interiors]
        simplified_poly = Polygon(outer, holes=inners)
        
        self.poly = orient(simplified_poly)

        self._mesh_v, self._mesh_f = triangulate_polygon(
            self.poly, engine='earcut')

        self.shell_bisector = self.compute_bisector(self.poly.exterior.coords)
        self.holes_bisector = [self.compute_bisector(
            hole.coords) for hole in self.poly.interiors]

    def compute_bisector(self, poly_v: CoordinateSequence, holes: bool = False):
        vert_i = np.array(poly_v)
        vec_a = vert_i[1:] - vert_i[0:-1]
        vec_b = -1*np.roll(vec_a, shift=1, axis=0)
        direction = np.cross(vec_a, vec_b, axis=1)

        # assume no colinear segment pair
        dir_sign = direction >= 0

        # flip direction for angle bigger than 180
        sign = np.where(dir_sign, np.ones(
            direction.shape[0]),  -1*np.ones(direction.shape[0]))[:, np.newaxis]
        bisector = LA.norm(vec_b, axis=1)[
            :, np.newaxis]*vec_a + LA.norm(vec_a, axis=1)[:, np.newaxis]*vec_b

        return sign*bisector/LA.norm(bisector, axis=1)[:, np.newaxis]

    def show_velocity_outer(self):
        polygon = load_path(self.poly)
        vert_1 = np.array(self.poly.exterior.coords)[:-1]
        vert_2 = vert_1 + self.shell_bisector
        lines = np.nan_to_num(np.stack((vert_1, vert_2), axis=1))
        shapely_lines = load_path(lines)

        (shapely_lines + polygon).show()

    def show_velocity_inner(self, idx: int = 0):
        polygon = load_path(self.poly)
        vert_1 = np.array(self.poly.interiors[idx].coords)[:-1]
        vert_2 = vert_1 + self.holes_bisector[idx]
        lines = np.nan_to_num(np.stack((vert_1, vert_2), axis=1))
        shapely_lines = load_path(lines)
        (shapely_lines + polygon).show()

    def show_mesh(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(self._mesh_v[:, 0], self._mesh_v[:, 1], np.zeros(
            self._mesh_v.shape[0]), triangles=self._mesh_f)
        plt.show()
