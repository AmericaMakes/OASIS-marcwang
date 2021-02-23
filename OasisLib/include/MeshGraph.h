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
BOOST_GEOMETRY_REGISTER_POINT_3D(pt, double, bg::cs::cartesian, data()[0], data()[1], data()[2]);
using c_range = bg::model::box<pt>;
using value = std::pair<c_range, MeshHalfedges::Halfedge>;
using rtree = bgi::rtree<value, bgi::quadratic<16> >;

namespace OasisLib{

    class MeshHeightSlicer{

        public:
            size_t nb_nodes;
            MeshHeightSlicer(std::shared_ptr<Mesh> m_ptr);
            vector<value> get_layer(double z);
            //void clip_cell(vector<value> &target_facet, Mesh &m_out, double z);
            //edge_pair clip_cell_by_plane(index_t f, const Mesh *mesh, double z);

        protected:
            rtree mesh_tree;
            std::shared_ptr<Mesh> mesh_ptr;
        private:
            void initialize_rtree();
    };

}
#endif