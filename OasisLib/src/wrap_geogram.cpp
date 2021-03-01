#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_remesh.h>
#include <geogram/basic/attributes.h>

#include "GeogramBase.h"
#include "GeogramVoronoi.h"
#include "GeogramCaster.h"
#include "MeshGraph.h"

using namespace OasisLib;
using namespace GEO;
namespace py = pybind11;

PYBIND11_MODULE(OasisLib, m)
{
    m.doc() = "python wrap geogram mesh";
    GeogramBase init;

    py::enum_<MeshElementsFlags>(m, "MeshElementsFlags", py::arithmetic())
        .value("MESH_NONE", MeshElementsFlags::MESH_NONE)
        .value("MESH_VERTICES", MeshElementsFlags::MESH_VERTICES)
        .value("MESH_FACETS", MeshElementsFlags::MESH_FACETS)
        .value("MESH_EDGES", MeshElementsFlags::MESH_EDGES)
        .value("MESH_CELLS", MeshElementsFlags::MESH_CELLS)
        .value("MESH_ALL_ELEMENTS", MeshElementsFlags::MESH_ALL_ELEMENTS)
        .value("MESH_FACET_CORNERS", MeshElementsFlags::MESH_FACET_CORNERS)
        .value("MESH_CELL_CORNERS", MeshElementsFlags::MESH_CELL_CORNERS)
        .value("MESH_CELL_FACETS", MeshElementsFlags::MESH_CELL_FACETS)
        .value("MESH_ALL_SUBELEMENTS", MeshElementsFlags::MESH_ALL_SUBELEMENTS)
        .export_values();

    py::class_<Mesh, std::shared_ptr<Mesh>>(m, "GeoMesh")
        .def(py::init<geo_index_t, bool>(),
             py::arg("dimension") = 3, py::arg("single_precision") = false)
        .def("clear", &Mesh::clear,
             py::arg("keep_attributes") = true, py::arg("keep_memory") = false)
        .def("show_stats", &Mesh::show_stats,
             py::arg("tag") = "Mesh")
        .def("copy", &Mesh::copy,
             py::arg("Mesh"), py::arg("copy_attributes") = true,
             py::arg("what") = MeshElementsFlags::MESH_ALL_ELEMENTS)
        .def("assert_is_valid", &Mesh::assert_is_valid)
        .def("get_attributes", &Mesh::get_attributes)
        .def("get_scalar_attributes", &Mesh::get_scalar_attributes)
        .def("get_vector_attributes", &Mesh::get_vector_attributes,
             py::arg("max_dim") = 0)
        .def("nb_subelements_types", &Mesh::nb_subelements_types)
        .def_readonly("cells", &Mesh::cells)
        .def_readonly("vertices", &Mesh::vertices)
        .def_readonly("facets", &Mesh::facets)
        .def_readonly("facet_corners", &Mesh::facet_corners)
        .def("get_facet_vertex", [](const Mesh &m_in, index_t f)
        {
            auto c_begin = m_in.facets.corners_begin(f);
            auto c_end = m_in.facets.corners_end(f);
            std::vector<index_t> v_id;
            for(auto i = c_begin; i != c_end; ++i)
            {
                auto v = m_in.facet_corners.vertex(i);
                v_id.push_back(v);
            }
            return v_id;
        });
    
    py::class_<MeshCells>(m, "MeshCells")
        .def("compute_borders", py::overload_cast<>(&MeshCells::compute_borders));
    
    py::class_<MeshVertices>(m, "MeshVertices")
        .def("point", py::overload_cast<index_t>(&MeshVertices::point, py::const_ ))
        .def("nb", &MeshVertices::nb);
    
    py::class_<MeshFacetCornersStore>(m, "MeshFacetCornersStore")
        .def("nb", &MeshVertices::nb)
        .def("vertex", &MeshFacetCornersStore::vertex)
        .def("adjacent_facet", &MeshFacetCornersStore::adjacent_facet);

    py::class_<MeshFacets>(m, "MeshFacets")
        .def("nb", &MeshFacets::nb)
        .def("get_cell_id", [](const MeshFacets &f){
            Attribute<int> cell_id(f.attributes(), "cell_id");
            auto size = cell_id.size();
            std::vector<int> a_cell;
            for(auto i = 0 ; i < size ; ++i)
            {
                a_cell.push_back(cell_id[i]);   
            }
            return a_cell;
        })
        .def("adjacent", &MeshFacets::adjacent)
        .def("vertex", &MeshFacets::vertex)
        .def("corners_begin", &MeshFacets::corners_begin)
        .def("corners_end", &MeshFacets::corners_end)
        .def("nb_corners", &MeshFacets::nb_corners)
        .def("nb_vertices", &MeshFacets::nb_vertices)
        .def("connect", &MeshFacets::connect)
        .def("nb_vertices", &MeshFacets::nb_vertices)
        .def("next_corner_around_facet", &MeshFacets::next_corner_around_facet);    

    py::enum_<MeshAttributesFlags>(m, "MeshAttributesFlags", py::arithmetic())
        .value("MESH_NO_ATTRIBUTES", MeshAttributesFlags::MESH_NO_ATTRIBUTES)
        .value("MESH_VERTEX_REGION", MeshAttributesFlags::MESH_VERTEX_REGION)
        .value("MESH_VERTEX_TEX_COORD", MeshAttributesFlags::MESH_VERTEX_TEX_COORD)
        .value("MESH_VERTEX_COLOR", MeshAttributesFlags::MESH_VERTEX_COLOR)
        .value("MESH_FACET_REGION", MeshAttributesFlags::MESH_FACET_REGION)
        .value("MESH_CELL_REGION", MeshAttributesFlags::MESH_CELL_REGION)
        .value("MESH_ALL_ATTRIBUTES", MeshAttributesFlags::MESH_ALL_ATTRIBUTES)
        .export_values();

    py::class_<MeshIOFlags>(m, "MeshIOFlags")
        .def(py::init<>())
        .def("dimension", &MeshIOFlags::dimension)
        .def("attributes", &MeshIOFlags::attributes)
        .def("set_attributes", &MeshIOFlags::set_attributes)
        .def("set_attribute", &MeshIOFlags::set_attribute)
        .def("has_attribute", &MeshIOFlags::has_attribute)
        .def("elements", &MeshIOFlags::elements)
        .def("set_elements", &MeshIOFlags::set_elements)
        .def("set_element", &MeshIOFlags::set_element)
        .def("has_element", &MeshIOFlags::has_element);

    m.def("mesh_load",
          py::overload_cast<const std::string &, Mesh &, const MeshIOFlags &>(&mesh_load),
          "load mesh", py::arg("filename"),
          py::arg("M"), py::arg("ioflags") = MeshIOFlags());
    
    m.def("mesh_save", 
        py::overload_cast<const Mesh &, const std::string &, const MeshIOFlags &>(&mesh_save), 
        "save mesh", py::arg("M"), py::arg("filename"), 
        py::arg("ioflags") = MeshIOFlags());
    
    py::enum_<MeshRepairMode>(m, "MeshRepairMode", py::arithmetic())
        .value("MESH_REPAIR_TOPOLOGY", MeshRepairMode::MESH_REPAIR_TOPOLOGY)
        .value("MESH_REPAIR_COLOCATE", MeshRepairMode::MESH_REPAIR_COLOCATE)
        .value("MESH_REPAIR_DUP_F", MeshRepairMode::MESH_REPAIR_DUP_F)
        .value("MESH_REPAIR_RECONSTRUCT", MeshRepairMode::MESH_REPAIR_RECONSTRUCT)
        .value("MESH_REPAIR_QUIET", MeshRepairMode::MESH_REPAIR_QUIET)
        .value("MESH_REPAIR_DEFAULT", MeshRepairMode::MESH_REPAIR_DEFAULT)
        .export_values();
    
    m.def("mesh_repair", &mesh_repair, 
        py::arg("M"), py::arg("mode") = MeshRepairMode::MESH_REPAIR_DEFAULT,
        py::arg("colocate_epsilon") = 0.0);
    
    m.def("mesh_tetrahedralize", &mesh_tetrahedralize, 
        py::arg("M"), py::arg("preprocess") = true, py::arg("refine") = false, 
        py::arg("quality") = 2.0, py::arg("keep_regions") = false);
    
    m.def("polyhedral_mesher", &polyhedral_mesher,
        py::arg("M_in"), py::arg("M_out"), 
        py::arg("nb_points") = 1000, py::arg("simplify") = "tets_voronoi_boundary",
        py::arg("angle_threshold") = 0.001, py::arg("nd_iter_lloyd") = 5,
        py::arg("nb_iter_newton") = 30, py::arg("tessallate_non_convex") = false,
        py::arg("poly_cell_shrinks") = 0.0);
    
    m.def("remesh_smooth", &remesh_smooth, 
        py::arg("M_in"), py::arg("M_out"), 
        py::arg("nb_points") = 10000, py::arg("dim") = 0,
        py::arg("nb_Lloyd_iter") = 5, py::arg("nb_Newton_iter") = 30,
        py::arg("Newton_m") = 7);

    py::class_<MeshHeightSlicer>(m, "MeshHeightSlicer")
        .def(py::init<std::shared_ptr<Mesh>>())
        .def("get_layer", &MeshHeightSlicer::get_layer);
}
