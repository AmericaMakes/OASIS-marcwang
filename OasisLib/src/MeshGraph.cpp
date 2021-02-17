#include "MeshGraph.h"
#include <geogram/basic/common.h>

#include <vector>
#include <limits>

namespace OasisLib
{
    rtree make_rtree(const Mesh &m_in)
    {
        rtree tree_inst;
        index_t nb_cell = m_in.facets.nb();

        
        for (index_t f = 0; f < m_in.facets.nb(); f++)
        {
            float z_min = std::numeric_limits<float>::infinity();
            float z_max = -std::numeric_limits<float>::infinity();

            for (index_t c = m_in.facets.corners_begin(f);
                 c < m_in.facets.corners_end(f); c++)
            {
                index_t g_index = m_in.facets.vertex(f, c);
                vec3 v_pos = m_in.vertices.point(g_index);

                if (v_pos[2] > z_max){
                    z_max = v_pos[2];
                }

                if (v_pos[2] < z_min){
                    z_min = v_pos[2];
                }

                c_range bx(pt(z_min), pt(z_max));
                //tree_inst.insert(std::make_pair(bx, f));

            }
        }
        
        return tree_inst;
    }

}