
from shapely.geometry.polygon import Polygon, orient
from trimesh.path.simplify import merge_colinear

def clean_polygon(poly : Polygon , precision : float = 1e-3):
        # remove colinear segment
    outer = merge_colinear(poly.exterior.coords, scale=precision)
    inners = [merge_colinear(hole.coords, scale=precision)
                for hole in poly.interiors]
    simplified_poly = Polygon(outer, holes=inners)

    return orient(simplified_poly)