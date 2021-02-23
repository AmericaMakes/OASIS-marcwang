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
        MeshHalfedges mhl(*this->mesh_ptr);
        
        index_t nb_f = this->mesh_ptr->facets.nb();
        for (index_t f_i = 0; f_i < nb_f; ++f_i)
        {
            for (auto c_i = this->mesh_ptr->facets.corners_begin(f_i);
                c_i < this->mesh_ptr->facets.corners_end(f_i); ++c_i)
            {
                MeshHalfedges::Halfedge hl(f_i, c_i);

                auto p_1_id = this -> mesh_ptr->facet_corners.vertex(c_i);
                auto c_i_1 = this -> mesh_ptr->facets.next_corner_around_facet(f_i, c_i);
                
                auto p_2_id = this -> mesh_ptr->facet_corners.vertex(c_i_1);

                auto p_1 = this -> mesh_ptr->vertices.point(p_1_id);
                auto p_2 = this -> mesh_ptr->vertices.point(p_2_id);
                
                vecng<3, double> min_p;
                vecng<3, double> max_p;

                for(int i = 0; i < 3; ++i)
                {
                    if(p_1[i] >= p_2[i])
                    {
                        min_p[i] = p_2[i];
                        max_p[i] = p_1[i];
                    }else
                    {
                        min_p[i] = p_1[i];
                        max_p[i] = p_2[i];
                    }
                }

                c_range bx(min_p,max_p);
                this->mesh_tree.insert(std::make_pair(bx, hl));
            }
        }
    }

    vector<value> MeshHeightSlicer::get_layer(double z)
    {
        vector<value> results;
        auto bds = this->mesh_tree.bounds();
        auto min_v = bds.min_corner();
        auto max_v = bds.max_corner();

        pt l_pt({min_v.get<0>(), min_v.get<1>(), z});
        pt m_pt({max_v.get<0>(), max_v.get<1>(), z});

        c_range bbx(l_pt, m_pt);

        this->mesh_tree.query(bgi::intersects(bbx), std::back_inserter(results));
        //Mesh m_out;
        //clip_cell(results, m_out, z);

        return results;
    }
    /*
    void MeshHeightSlicer::clip_cell(vector<value> &target_facet, Mesh &m_out, double z)
    {   
        auto nb_e = target_facet.size();
        //auto e_idx = m_out.create_edges(nb_e);
        std::vector<edge_pair> e_v;
        for (auto it = target_facet.cbegin(); it != target_facet.cend(); ++it)
        {
            auto pair = *it;
            auto c_i = pair.second;
            //Polygon result;
            //auto e_pair = clip_poly_by_plane(c_i, this->mesh_ptr.get(), z);
            //e_v.push_back(e_pair);

            //auto p1 = e_pair.first;
            //auto p2 = e_pair.second;

            //auto p1_idx = m_out.vertices.create_vertex(p1.data());
            //auto p2_idx = m_out.vertices.create_vertex(p2.data());

            //m_out.edges.set_vertex(e_idx, 0, p1_idx);
            //m_out.edges.set_vertex(e_idx, 1, p2_idx);
            //e_idx++;
        }

    }
    
    index_t MeshHeightSlicer::clip_cell_by_plane(index_t c_i, const Mesh *mesh, Mesh &out,double z)
    {
        vecng<3, double> p0{0.0, 0.0, z};
        vecng<3, double> n{0.0, 0.0, 1.0};

        auto mhl = MeshHalfedges(*mesh);
        auto nb_facets_out = out.facets.nb();
        auto nb_f = mesh -> cells.nb_facets(c_i);

        vector<index_t> poly;
        for(auto f = 0; k < nb_f; ++f)
        {
            auto nb_v = mesh -> cells.facet_nb_vertices(c_i, f);

            Halfedge hl(f, 0);

            int count = 0;
            auto c_idx = mesh -> cells.facet_corner(c_i, hl.facet, hl.corner);
            index_t current_v_idx = mesh->facet_corners.vertex(c_idx);
            auto c_pos = mesh->vertices.point(current_v_idx);

            mhl.move_to_next_around_facet(hl);          
            auto n_idx = mesh -> cells.facet_corner(c_i, hl.facet, hl.corner);
            index_t n_v_idx = mesh->facet_corners.vertex(n_idx);
            auto n_pos = mesh->vertices.point(n_v_idx);

            auto is_intersect = PCK::dot_compare_3d(n.data(), n_pos.data(), c_pos.data());
            auto is_parrallel = PCK::dot_compare_3d(n.data(), p0.data(), c_pos.data());

            // if sign not identical then intersection
            if (is_intersect != 0)
            {
                auto inter_p = c_pos + dot((p0 - c_pos), n)/dot((n_pos - c_pos),n)*(n_pos - c_pos);
                auto intersect_v_idx = mesh -> create_vertex(inter_p.data());
                poly.push_back(intersect_v_idx);
            }else{
                if (is_parrallel == 0 && c == 0)
                {
                    
                }
            }
        }
        return 0;
    }
    */
}