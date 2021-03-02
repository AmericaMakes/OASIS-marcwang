import pyclipper
from shapely.geometry import Polygon
from typing import List

from oasis.util import clean_polygon, convert_to_clipper, convert_from_clipper

class ClipperOffsetting():

    def __init__(self, poly : Polygon, clipper_scale : int = 1000) -> None:
        self.poly = clean_polygon(poly)
        self.offsetter = pyclipper.PyclipperOffset()
        self.clipper_scale = clipper_scale
        cl_path = convert_to_clipper(self.poly, self.clipper_scale)
        self.offsetter.AddPaths(cl_path, pyclipper.JT_ROUND, pyclipper.ET_CLOSEDPOLYGON)
    
    def get_offset(self, dist : float, n_contour: int = 1) -> List[Polygon]:
        result = []
        for i in range(n_contour):
            c_path = self.offsetter.Execute(-1*dist*i*self.clipper_scale)
            poly = convert_from_clipper(c_path, self.clipper_scale)
            result.append(poly)
        return result