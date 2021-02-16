#include "GeogramVoronoi.h"

#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/RVD_mesh_builder.h>

namespace OasisLib{

int polyhedral_mesher(Mesh &M_in, Mesh &M_out,
                      int nb_points,
                      const string& simplify,
                      double angle_threshold,
                      index_t nb_iter_lloyd,
                      index_t nb_iter_newton,
                      bool tessallate_non_convex,
                      float poly_cell_shrinks,
                      bool generate_ids)
{
    Logger::div("Polyhedral meshing");

    if (M_in.cells.nb() == 0)
    {
        Logger::out("Poly") << "Mesh is not a volume" << std::endl;
        Logger::out("Poly") << "Trying to tetrahedralize" << std::endl;
        if (!mesh_tetrahedralize(M_in))
        {
            return 0;
        }
        M_in.cells.compute_borders();
    }

    index_t dim = M_in.vertices.dimension();

    CentroidalVoronoiTesselation CVT(&M_in, coord_index_t(dim));
    CVT.set_volumetric(true);
    Logger::div("Generate random samples");

    CVT.compute_initial_sampling(nb_points);

    Logger::div("Optimize sampling");
    {
        ProgressTask progress("Lloyd", nb_iter_lloyd);
        CVT.set_progress_logger(&progress);
        CVT.Lloyd_iterations(nb_iter_lloyd);
    }
    {
        ProgressTask progress("Newton", nb_iter_newton);
        CVT.set_progress_logger(&progress);
        CVT.Newton_iterations(nb_iter_newton);
    }
    CVT.set_progress_logger(nullptr);

    CVT.RVD()->set_exact_predicates(true);
    {
        BuildRVDMesh callback(M_out);
        if (simplify == "tets_voronoi_boundary")
        {
            callback.set_simplify_boundary_facets(true, angle_threshold);
        }
        else if (simplify == "tets_voronoi")
        {
            callback.set_simplify_voronoi_facets(true);
        }
        else if (simplify == "tets")
        {
            callback.set_simplify_internal_tet_facets(true);
        }
        else if (simplify == "none")
        {
            callback.set_simplify_internal_tet_facets(false);
        }
        else
        {
            Logger::err("Poly")
                << simplify << " invalid cells simplification mode"
                << std::endl;
            return 0;
        }

        callback.set_tessellate_non_convex_facets(tessallate_non_convex);
        callback.set_shrink(poly_cell_shrinks);
        callback.set_generate_ids(generate_ids);
        CVT.RVD()->for_each_polyhedron(callback);
    }

    return 1;
}

}