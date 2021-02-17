#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_remesh.h>
#include <geogram/mesh/mesh_tetrahedralize.h>

#include <boost/geometry/index/rtree.hpp>

#include "MeshGraph.h"

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/common.h>

using namespace GEO;

TEST_CASE("Test slice query", "[Query]"){
    std::string path("./test_artifact/3DBenchy.stl");

    initialize();
    CmdLine::import_arg_group("standard");
    CmdLine::import_arg_group("pre");
    CmdLine::import_arg_group("remesh");
    CmdLine::import_arg_group("algo");
    CmdLine::import_arg_group("post");
    CmdLine::import_arg_group("opt");
    CmdLine::import_arg_group("co3ne");
    CmdLine::import_arg_group("tet");
    CmdLine::import_arg_group("poly");

    Mesh m_in;
    Mesh m_smooth;
    REQUIRE(mesh_load(path, m_in));
    mesh_repair(m_in);

    index_t nb_points = 1000;
    remesh_smooth(m_in, m_smooth, nb_points);

    REQUIRE(mesh_tetrahedralize(m_smooth));

    SECTION("test rtree construction"){
        OasisLib::make_rtree(m_smooth);
    }

}