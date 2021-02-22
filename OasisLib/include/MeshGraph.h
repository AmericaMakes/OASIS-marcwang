#ifndef MESHGRAPH_H
#define MESHGRAPH_H

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>
#include <geogram/basic/numeric.h>

using namespace GEO;
using namespace GEOGen;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using pt = bg::model::point<float, 3, bg::cs::cartesian>;
using c_range = bg::model::box<pt>;
using value = std::pair<c_range, index_t>;
using rtree = bgi::rtree<value, bgi::quadratic<16> >;
using cell_id_map = std::map<index_t, std::vector<index_t> >;

namespace OasisLib{

    class MeshHeightSlicer{

        public:
            size_t nb_nodes;
            MeshHeightSlicer(std::shared_ptr<Mesh> m_ptr);
            vector<value> get_layer(double z);
            void clip_cell(vector<value> &target_cell, Mesh &m_out, cell_id_map &id_map, double z);
            void clip_poly_by_plane(Polygon &target, Polygon& result, const Mesh *mesh, double *pi, double *pj);

        protected:
            rtree mesh_tree;
            std::shared_ptr<Mesh> mesh_ptr;
        private:
            void initialize_rtree();
    };

}
#endif