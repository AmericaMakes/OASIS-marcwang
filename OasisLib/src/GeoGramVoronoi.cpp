#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/process.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/geometry_nd.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_tetrahedralize.h>

#include <geogram/delaunay/LFS.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/RVD_mesh_builder.h>

#include <typeinfo>
#include <algorithm>

using namespace GEO;

int polyhedral_mesher(
    const std::string &input_filename, std::string output_filename)
{
    Mesh M_in;
    Mesh M_out;
    Mesh M_points;

    Logger::div("Polyhedral meshing");

    if (!mesh_load(input_filename, M_in))
    {
        return 1;
    }

    if (M_in.cells.nb() == 0)
    {
        Logger::out("Poly") << "Mesh is not a volume" << std::endl;
        Logger::out("Poly") << "Trying to tetrahedralize" << std::endl;
        if (!mesh_tetrahedralize(M_in))
        {
            return 1;
        }
        M_in.cells.compute_borders();
    }

    index_t dim = M_in.vertices.dimension();
    index_t spec_dim = CmdLine::get_arg_uint("poly:embedding_dim");
    if (spec_dim != 0 && spec_dim <= dim)
    {
        dim = spec_dim;
    }

    CentroidalVoronoiTesselation CVT(&M_in, coord_index_t(dim));
    CVT.set_volumetric(true);

    if (CmdLine::get_arg("poly:points_file") == "")
    {

        Logger::div("Generate random samples");

        CVT.compute_initial_sampling(
            CmdLine::get_arg_uint("remesh:nb_pts"));

        Logger::div("Optimize sampling");

        try
        {
            index_t nb_iter = CmdLine::get_arg_uint("opt:nb_Lloyd_iter");
            ProgressTask progress("Lloyd", nb_iter);
            CVT.set_progress_logger(&progress);
            CVT.Lloyd_iterations(nb_iter);
        }
        catch (const TaskCanceled &)
        {
        }

        try
        {
            index_t nb_iter = CmdLine::get_arg_uint("opt:nb_Newton_iter");
            ProgressTask progress("Newton", nb_iter);
            CVT.set_progress_logger(&progress);
            CVT.Newton_iterations(nb_iter);
        }
        catch (const TaskCanceled &)
        {
        }

        CVT.set_progress_logger(nullptr);
    }
    else
    {
        if (!mesh_load(CmdLine::get_arg("poly:points_file"), M_points))
        {
            return 1;
        }
        CVT.delaunay()->set_vertices(
            M_points.vertices.nb(), M_points.vertices.point_ptr(0));
    }

    CVT.RVD()->set_exact_predicates(true);
    {
        BuildRVDMesh callback(M_out);
        std::string simplify = CmdLine::get_arg("poly:simplify");
        if (simplify == "tets_voronoi_boundary")
        {
            double angle_threshold =
                CmdLine::get_arg_double("poly:normal_angle_threshold");
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
        }
    
    callback.set_tessellate_non_convex_facets(
        CmdLine::get_arg_bool("poly:tessellate_non_convex_facets"));
    callback.set_shrink(CmdLine::get_arg_double("poly:cells_shrink"));
    callback.set_generate_ids(
        CmdLine::get_arg_bool("poly:generate_ids") ||
        FileSystem::extension(output_filename) == "ovm");
    CVT.RVD()->for_each_polyhedron(callback);
    }
}