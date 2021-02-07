from shapely.geometry import Polygon, MultiLineString, LineString
from shapely import affinity
import trimesh
from trimesh.voxel.creation import voxelize
from trimesh.intersections import plane_lines
from enum import IntEnum
import numpy as np

from scipy import stats
from scipy.spatial import Voronoi, ConvexHull

import networkx as nx
from rtree import index


class HatchPattern(IntEnum):
    Straight = 0


class PointSampler():
    def __init__(self, mesh: trimesh.Trimesh, size: float) -> None:
        self._vox = voxelize(mesh, size, method='ray')
        self.n_box = len(self._vox.points)
        self._size = size

    def generate_point(self, n_point: int):
        grid_idx = np.random.randint(0, self.n_box, size=n_point)
        mid_point = np.array(self._vox.points[grid_idx])
        offset = stats.uniform.rvs(-0.5 * self._size,
                                   self._size, size=n_point*3).reshape(n_point, 3)
        return mid_point + offset


class VoronoiInfill():
    def __init__(self, mesh: trimesh.Trimesh, r_max: float, n_cell: int = 1000) -> None:
        self._mesh = mesh
        self._sampler = PointSampler(mesh, r_max)

        pt = self._sampler.generate_point(n_cell)
        self.vor = Voronoi(pt)
        self._vor_tree = self._build_rtree(self.vor)
        #self.vor_cell = self.get_voronoi_cell(self.vor)

    def _build_rtree(self, vor):
        p = index.Property()
        p.dimension = 3
        r_tree = index.Index(properties=p)
        for i, simplex in enumerate(vor.regions):
            if len(simplex) > 0:
                simplex = np.array(simplex)
                mask = simplex >= 0
                r_v = vor.vertices[simplex[mask]]
                coord = (np.min(r_v[:, 0]), np.min(r_v[:, 1]), np.min(r_v[:, 2]),
                        np.max(r_v[:, 0]), np.max(r_v[:, 1]), np.max(r_v[:, 2]))
                r_tree.insert(i, coord)
        return r_tree

    def edge_plane_intersect(self, z, edges):
        i = 0
        proj_line = np.zeros((2, 3))
        for e in edges:
            i_point, valid = plane_lines((0, 0, z), (0, 0, 1), e)
            if valid:
                proj_line[i] = i_point
                i += 1
            if i == 1:
                break
        return proj_line

    def get_slice(self, z_h: float):
        vor_g = nx.Graph()
        vert_g = nx.Graph()
        bounds = self._mesh.bounding_box
        #hits = self._vor_tree.intersection()

def straight_hatch(contour: Polygon, offset: float, angle: float):
    center_rot = contour.centroid
    r_contour = affinity.rotate(contour, angle, origin=center_rot)
    bounds = r_contour.bounds
    n_contour = np.ceil((bounds[2] - bounds[0])/offset).astype(int)
    ls_line = [((i, bounds[1]), (i, bounds[3]))
               for i in np.linspace(bounds[0], bounds[2], n_contour)]
    hatch = MultiLineString(ls_line)
    clipped_hatch = hatch.intersection(r_contour)
    filtered = MultiLineString([geo for geo in list(
        clipped_hatch) if isinstance(geo, LineString)])

    return affinity.rotate(filtered, -angle, origin=center_rot)
