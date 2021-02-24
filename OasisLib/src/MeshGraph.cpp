#include "MeshGraph.h"
#include <geogram/basic/common.h>
#include <geogram/numerics/predicates.h>
#include <limits>

namespace OasisLib
{
    MeshHeightSlicer::MeshHeightSlicer(std::shared_ptr<Mesh> m_ptr)
    {
        this->mesh_ptr = m_ptr;
        this->initialize_rtree();
        this->nb_nodes = this->mesh_tree.size();
    }

    void MeshHeightSlicer::initialize_rtree()
    {
        index_t nb_c = this->mesh_ptr->cells.nb();
        for (index_t c = 0; c < nb_c; ++c)
        {
            index_t nb_e = this -> mesh_ptr->cells.nb_edges(c);

            for(auto e_i = 0; e_i < nb_e; ++e_i){
                
                auto e0_v = this->mesh_ptr->cells.edge_vertex(c, e_i, 0);
                auto p_1 = this->mesh_ptr->vertices.point(e0_v);

                auto e1_v = this->mesh_ptr->cells.edge_vertex(c, e_i, 1);
                auto p_2 = this->mesh_ptr->vertices.point(e1_v);

                vecng<3, double> min_p;
                vecng<3, double> max_p;

                for (int i = 0; i < 3; ++i)
                {
                    if (p_1[i] >= p_2[i])
                    {
                        min_p[i] = p_2[i];
                        max_p[i] = p_1[i];
                    }
                    else
                    {
                        min_p[i] = p_1[i];
                        max_p[i] = p_2[i];
                    }
                }

                c_range bx(min_p, max_p);
                fc_p ids{c, e_i};
                this->mesh_tree.insert(std::make_pair(bx, ids));
                
            }   
            
        }
    }

    std::shared_ptr<Mesh> MeshHeightSlicer::get_layer(double z)
    {
        vector<value> results;
        auto bds = this->mesh_tree.bounds();
        auto min_v = bds.min_corner();
        auto max_v = bds.max_corner();

        pt l_pt({min_v.get<0>(), min_v.get<1>(), z});
        pt m_pt({max_v.get<0>(), max_v.get<1>(), z});

        c_range bbx(l_pt, m_pt);

        this->mesh_tree.query(bgi::intersects(bbx), std::back_inserter(results));
        
        auto m_out = std::make_shared<Mesh>();
        clip_cell(results, *m_out, z);

        return m_out;
    }

    void MeshHeightSlicer::clip_cell(vector<value> &target_facet, Mesh &m_out, double z)
    {
        std::map<index_t, std::map<index_t, std::set<index_t>>> processed_c;

        for (auto it = target_facet.cbegin(); it != target_facet.cend(); ++it)
        {
            auto ids = it -> second;

            auto c_it = processed_c.find(ids.first);

            if(c_it == processed_c.cend())
            {
                std::map<index_t, std::set<index_t>> poly_e;
                processed_c.insert({ids.first, poly_e});
            }
            

            auto f = this->mesh_ptr->cells.edge_adjacent_facet(ids.first, ids.second, 0);
            auto f_opp = this->mesh_ptr->cells.edge_adjacent_facet(ids.first, ids.second, 1);

            auto p_1_id = this->mesh_ptr->cells.edge_vertex(ids.first, ids.second, 0);
            auto p_1 = this->mesh_ptr->vertices.point(p_1_id);

            auto c_1_sign = geo_sgn(p_1[2] - z);

            auto p_2_id = this->mesh_ptr->cells.edge_vertex(ids.first, ids.second, 1);
            auto p_2 = this->mesh_ptr->vertices.point(p_2_id);

            auto c_2_sign = geo_sgn(p_2[2] - z);
            index_t inter_id;
            bool intersection = false;
            if (c_1_sign != c_2_sign && c_1_sign != 0)
            {
                double d = (z - p_1[2]) / (p_2[2] - p_1[2]);
                auto inter_p = p_1 + d * (p_2 - p_1);
                inter_id = m_out.vertices.create_vertex(inter_p.data());
                intersection = true;
            }

            if (c_1_sign == 0)
            {
                inter_id = m_out.vertices.create_vertex(p_1.data());
                intersection = true;
            } 

            if(intersection)
            {
                auto poly = processed_c[ids.first];

                auto f_it = poly.find(f);

                if(f_it == poly.cend())
                {
                    std::set<index_t> edge;
                    poly.insert({f, edge});
                }

                poly[f].insert(inter_id);

                auto f_opp_it = poly.find(f_opp);

                if(f_opp_it == poly.cend())
                {
                    std::set<index_t> edge;
                    poly.insert({f_opp, edge});
                }

                poly[f_opp].insert(inter_id);

            }
            
        }

        std::vector<index_t> cell_map;
        for(auto it = processed_c.cbegin(); it != processed_c.cend(); ++it)
        {
            cell_map.push_back(it -> first);
            std::vector<index_t> vertex_list;
            auto edge_list = it -> second;

            for(auto it_e = edge_list.cbegin(); it_e != edge_list.cend(); ++it_e)
            {
                auto edge = it_e -> second;
                if(edge.size() > 1 && edge.size() < 2)
                {
                    auto e_i = m_out.edges.create_edge();
                    int count = 0;
                    for(auto e_it = edge.cbegin(); e_it != edge.cend(); ++e_it)
                    {
                        auto v_id = *e_it;
                        m_out.edges.set_vertex(e_i, count, v_id);
                        count++;
                    }

                }
            }
            
        }
        
    }
}