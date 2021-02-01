import pyclipper
from shapely.geometry import Polygon
from shapely import affinity

from enum import IntEnum
import numpy as np

class HatchPattern(IntEnum):
    Straight = 0 

class Hatching():
    def __init__(self, style : int, clipper_scale : float = 1000) -> None:
        self.style = style
        self._scale = clipper_scale

    def straight_hatch(self, contour : Polygon, offset : float, angle : float):
        r_contour = affinity.rotate(contour, angle)
        bounds = r_contour.bounds
        n_contour = np.ceil((bounds[1] - bounds[0])/offset)
        hatch = [[(i,bounds[2]),(i,bounds[3])] for i in np.linspace(bounds[0], bounds[1], n_contour)]
        pyclipper.scale_to_clipper()
    
    
