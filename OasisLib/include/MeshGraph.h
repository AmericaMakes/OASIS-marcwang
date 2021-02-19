#ifndef MESHGRAPH_H
#define MESHGRAPH_H

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>

using namespace GEO;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using pt = bg::model::point<float, 3, bg::cs::cartesian>;
using c_range = bg::model::box<pt>;
using value = std::pair<c_range, index_t>;

using rtree = bgi::rtree<value, bgi::quadratic<16> >;

namespace OasisLib{

    class MeshHeightSlicer{

        public:
            size_t nb_nodes;

            MeshHeightSlicer(std::shared_ptr<Mesh> m_ptr);
            void get_layer(double z);

        protected:
            rtree mesh_tree;
            std::shared_ptr<Mesh> mesh_ptr;
        private:
            void initialize_rtree();

    };

}
#endif