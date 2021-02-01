import pyclipper
from shapely.geometry import Polygon, MultiLineString, LineString
from shapely import affinity

from enum import IntEnum
import numpy as np

from oasis.util import convert_from_clipper, convert_to_clipper


class HatchPattern(IntEnum):
    Straight = 0


class Hatching():
    def __init__(self, style: int, clipper_scale: float = 1000) -> None:
        self.style = style
        self._scale = clipper_scale

    def straight_hatch(self, contour: Polygon, offset: float, angle: float):
        center_rot = contour.centroid
        r_contour = affinity.rotate(contour, angle, origin = center_rot)
        bounds = r_contour.bounds
        n_contour = np.ceil((bounds[2] - bounds[0])/offset).astype(int)
        ls_line = [((i, bounds[1]), (i, bounds[3]))
                 for i in np.linspace(bounds[0], bounds[2], n_contour)]
        hatch = MultiLineString(ls_line)
        clipped_hatch = hatch.intersection(r_contour)
        filtered = MultiLineString([geo for geo in list(clipped_hatch) if isinstance(geo, LineString)])

        return affinity.rotate(filtered, -angle, origin = center_rot)
