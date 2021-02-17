#include "MeshGraph.h"
#include <geogram/basic/common.h>

#include <vector>
#include <limits>

namespace OasisLib
{

    void make_rtree(Mesh &m_in)
    {
        index_t nb_cell = m_in.facets.nb();

        auto m_tree = std::make_shared<mesh_tree>();
        no_iterator b = m_in.cells.begin();
        no_iterator e = m_in.cells.end();

        std::vector<value> result_s;

        for (index_t f = 0; f < m_in.facets.nb(); f++)
        {
            double z_min = std::numeric_limits<double>::infinity();
            double z_max = -std::numeric_limits<double>::infinity();

            for (index_t c = m_in.facets.corners_begin(f);
                 c < m_in.facets.corners_end(f); c++)
            {
                index_t g_index = m_in.facets.vertex(f, c);
                vec3 v_pos = m_in.vertices.point(g_index);
            }
        }
        /*
        while(b != e)
        {
            range<no_iterator> v = m_in.cells.corners(b);
            double z_min = std::numeric_limits<double>::infinity();
            double z_max = -std::numeric_limits<double>::infinity();

            for (i = range.begin();i < range.end(); i++){
                vec3 vertex = vs.point(i);
                if (vertex[2] > z_max){
                    z_max = vertex[2];
                }

                if (vertex[2] < z_min){
                    z_min = vertex[2];
                }
            }

            box bx(pt(z_min), pt(z_max));
            m_tree -> insert(value(bx, b));
            b++;
        }

        return m_tree;
        */
    }

} // namespace OasisLib