#ifndef MESHGRAPH_H
#define MESHGRAPH_H

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>
#include <geogram/basic/numeric.h>

using namespace GEO;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using pt = bg::model::point<float, 3, bg::cs::cartesian>;
using c_range = bg::model::box<pt>;
using value = std::pair<c_range, index_t>;
using rtree = bgi::rtree<value, bgi::quadratic<16> >;
using cell_id_map = std::map<index_t, std::vector<index_t> >;
using edge_pair = std::pair<vecng<3, double>,vecng<3, double>>;

namespace OasisLib{

    class MeshHeightSlicer{

        public:
            size_t nb_nodes;
            MeshHeightSlicer(std::shared_ptr<Mesh> m_ptr);
            vector<value> get_layer(double z);
            void clip_cell(vector<value> &target_facet, Mesh &m_out, double z);
            edge_pair clip_cell_by_plane(index_t f, const Mesh *mesh, double z);

        protected:
            rtree mesh_tree;
            std::shared_ptr<Mesh> mesh_ptr;
        private:
            void initialize_rtree();
    };

}
#endif