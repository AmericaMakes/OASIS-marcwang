#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_remesh.h>
#include <geogram/mesh/mesh_tetrahedralize.h>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "MeshGraph.h"

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/common.h>

#include <iostream>

using namespace GEO;
namespace bgi = boost::geometry::index;

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
    auto m_in = std::make_shared<Mesh>();
    auto m_smooth = std::make_shared<Mesh>();
    REQUIRE(mesh_load(path, *m_in));
    mesh_repair(*m_in);

    index_t nb_points = 10000;
    remesh_smooth(*m_in, *m_smooth, nb_points);

    REQUIRE(mesh_tetrahedralize(*m_smooth));

    auto nb_cell = m_smooth->cells.nb();
    SECTION("test rtree construction"){
        auto t = OasisLib::MeshHeightSlicer(m_smooth);
        REQUIRE(t.nb_nodes > 0);
        auto intersect_layer = t.get_layer(20);
        REQUIRE(intersect_layer.size() > 0);
    }

}