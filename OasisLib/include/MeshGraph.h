#ifndef MESHGRAPH_H
#define MESHGRAPH_H

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>

using namespace GEO;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<float, 1, bg::cs::cartesian> pt;
typedef bg::model::box<pt> box;
typedef std::pair<box, index_t> value;

typedef bgi::rtree<value, bgi::quadratic<16> > mesh_tree;

namespace OasisLib{
    void make_rtree(Mesh &m_in);
}
#endif