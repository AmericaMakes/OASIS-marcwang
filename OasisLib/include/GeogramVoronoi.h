#ifndef GEOGRAMVORONOI_H
#define GEOGRAMVORONOI_H

#include <string> 
#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>

using namespace std;
using namespace GEO;

namespace OasisLib{
    int polyhedral_mesher(Mesh &M_in, Mesh &M_out,
                        int nb_points = 1000,
                        const string& simplify = "tets_voronoi_boundary",
                        double angle_threshold = 0.001,
                        index_t nb_iter_lloyd = 5,
                        index_t nb_iter_newton = 30,
                        bool tessallate_non_convex = false,
                        float poly_cell_shrinks = 0.0);
}
#endif