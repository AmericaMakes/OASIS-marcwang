#ifndef MESHGRAPH_H
#define MESHGRAPH_H

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/register/point.hpp>

#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_halfedges.h>
#include <geogram/basic/numeric.h>

using namespace GEO;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using pt = vecng<3, double>;
using pt2d = vecng<2, double>;
BOOST_GEOMETRY_REGISTER_POINT_3D(pt, double, bg::cs::cartesian, data()[0], data()[1], data()[2]);
BOOST_GEOMETRY_REGISTER_POINT_2D(pt2d, double, bg::cs::cartesian, data()[0], data()[1]);
using c_range = bg::model::box<pt>;
using mpt2d = bg::model::multi_point<pt2d>;
using contour2d = bg::model::ring<pt2d>;
using value = std::pair<c_range, std::pair<index_t, index_t>>;
using rtree = bgi::rtree<value, bgi::quadratic<16>>;

using sh_mesh_ptr = std::shared_ptr<Mesh>;
using id_map = std::map<index_t, index_t>;

namespace OasisLib
{

    class MeshHeightSlicer
    {
    public:

        size_t nb_nodes;
        MeshHeightSlicer(sh_mesh_ptr m_ptr);
        std::pair<sh_mesh_ptr, id_map> get_layer(double z);

    protected:
        rtree mesh_tree;
        sh_mesh_ptr mesh_ptr;

    private:
        id_map clip_cell(vector<value> &target_facet, Mesh &m_out, double z);
        void initialize_rtree();
    };

}
#endif