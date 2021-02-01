import pyclipper
from shapely.geometry import Polygon

from oasis.util import clean_polygon

class ClipperOffsetting():

    def __init__(self, poly : Polygon, clipper_scale : int = 1000):
        self.poly = clean_polygon(poly)
        self.offsetter = pyclipper.PyclipperOffset()
        self.clipper_scale = clipper_scale
        pyclipper.scale_to_clipper(clipper_scale)
        ls_hole = [pyclipper.scale_to_clipper(list(hole.coords), clipper_scale) for hole in self.poly.interiors]
        ls_path = [pyclipper.scale_to_clipper(list(self.poly.exterior.coords), clipper_scale)] + ls_hole
        
        self.offsetter.AddPaths(ls_path, pyclipper.JT_ROUND, pyclipper.ET_CLOSEDPOLYGON)
    
    def get_offset(self, dist : float, n_contour: int = 1):
        result = []
        for i in range(n_contour):
            path = self.offsetter.Execute(-1*dist*i*self.clipper_scale)
            scaled_p = pyclipper.scale_from_clipper(path, self.clipper_scale)
            poly = Polygon(shell = scaled_p[0], holes=scaled_p[1:])
            result.append(poly)
        return result