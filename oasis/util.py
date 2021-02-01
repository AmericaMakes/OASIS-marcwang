
from shapely.geometry.polygon import Polygon, orient
from trimesh.path.simplify import merge_colinear
import pyclipper as pcl
from typing import List, Tuple


def clean_polygon(poly: Polygon, precision: float = 1e-3):
    # remove colinear segment
    outer = merge_colinear(poly.exterior.coords, scale=precision)
    inners = [merge_colinear(hole.coords, scale=precision)
              for hole in poly.interiors]
    simplified_poly = Polygon(outer, holes=inners)

    return orient(simplified_poly)


def convert_to_clipper(poly: Polygon, scale: float = 1000):
    ls_path = [pcl.scale_to_clipper(list(hole.coords), scale)
               for hole in poly.interiors]
    if len(ls_path) > 0:
        ls_path = [pcl.scale_to_clipper(
            list(poly.exterior.coords), scale)] + ls_path
    return ls_path


def convert_from_clipper(cli_path: List[Tuple[float, float]], scale: float = 1000):
    scaled_p = pcl.scale_from_clipper(cli_path, scale)
    poly = Polygon(shell=scaled_p[0], holes=scaled_p[1:])
    return poly
