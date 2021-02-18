#include "MeshGraph.h"
#include <geogram/basic/common.h>

#include <limits>

namespace OasisLib
{
    rtree make_rtree(const Mesh &m_in)
    {
        rtree tree_inst;
        index_t nb_cell = m_in.facets.nb();
        
        for (index_t f = 0; f < m_in.facets.nb(); ++f)
        {
            std::array<double, 3> min_v;
            min_v.fill(std::numeric_limits<float>::infinity());

            std::array<double, 3> max_v;
            max_v.fill(-std::numeric_limits<float>::infinity());

            for (index_t c = 0;
                 c < m_in.facets.nb_vertices(f); ++c)
            {
                index_t g_index = m_in.facets.vertex(f, c);
                auto v_pos = m_in.vertices.point(g_index);

                for(index_t i = 0; i < 3 ; i++)
                {
                    if( v_pos[i] < min_v[i]){
                        min_v[i] = v_pos[i];
                    }

                    if(v_pos[i] > max_v[i])
                    {
                        max_v[i] = v_pos[i];
                    }
                }
            }

            c_range bx(pt(min_v[0], min_v[1], min_v[2]), 
                            pt(max_v[0], max_v[1], max_v[2]));
                
            tree_inst.insert(std::make_pair(bx, f));
        }
        
        return tree_inst;
    }

}