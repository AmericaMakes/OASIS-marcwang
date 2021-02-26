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
using value = std::pair<c_range, index_t>;
using rtree = bgi::rtree<value, bgi::quadratic<16> >;

namespace OasisLib{

    class MeshHeightSlicer{
        public:
            //enum class EdgeIntersect{ No, One, Two};
            size_t nb_nodes;
            MeshHeightSlicer(std::shared_ptr<Mesh> m_ptr);
            std::shared_ptr<Mesh> get_layer(double z);
            void clip_cell(vector<value> &target_facet, Mesh &m_out, double z);
            bool IsFacet(index_t c, index_t f, double z);
            int IsIntersection(index_t c_1, index_t c_2, pt &inter_p, double z);

        protected:
            rtree mesh_tree;
            std::shared_ptr<Mesh> mesh_ptr;
        private:
            void initialize_rtree();
    };

}
#endif