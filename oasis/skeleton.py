import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA

from shapely.geometry.polygon import Polygon, orient
from shapely.coords import CoordinateSequence
from shapely.geometry import MultiLineString
from trimesh.creation import triangulate_polygon
from trimesh.path.exchange.load import load_path


class Skeletonization():
    def __init__(self, polygon: Polygon) -> None:
        self.poly = orient(polygon)
        self._mesh_v, self._mesh_f = triangulate_polygon(
            polygon, engine='earcut')

        self.inner = self.compute_bisector(self.poly.exterior.coords)
        self.outer = [self.compute_bisector(
            hole.coords, False) for hole in self.poly.interiors]

    def compute_bisector(self, poly_v: CoordinateSequence, ccw: bool = True):
        vert_i = np.array(poly_v)
        vert_i_1 = np.roll(vert_i, shift=1, axis=0)
        vec_a = vert_i_1 - vert_i
        vec_b = -1*np.roll(vec_a, shift=1, axis=0)
        direction = np.cross(vec_a, vec_b, axis=1)

        if (ccw):
            dir_sign = direction <= 0
        else:
            dir_sign = direction > 0

        sign = np.where(dir_sign, np.ones(
            direction.shape[0]),  -1*np.ones(direction.shape[0]))[:, np.newaxis]
        bisector = LA.norm(vec_b, axis=1)[
            :, np.newaxis]*vec_a + LA.norm(vec_a, axis=1)[:, np.newaxis]*vec_b

        return sign*bisector/LA.norm(bisector, axis=1)[:, np.newaxis]
    
    def show_velocity(self):
        polygon = load_path(self.poly)
        vert_1 = np.array(self.poly.exterior.coords)
        vert_2 = vert_1 + self.inner
        
        shapely_line = MultiLineString(np.nan_to_num(np.stack((vert_1, vert_2), axis=2)).tolist())

        lines = load_path(shapely_line)

        (lines + polygon).show()

    def show_mesh(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(self._mesh_v[:, 0], self._mesh_v[:, 1], np.zeros(
            self._mesh_v.shape[0]), triangles=self._mesh_f)
        plt.show()
