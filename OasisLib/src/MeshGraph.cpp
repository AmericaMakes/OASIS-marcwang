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
        index_t nb_cell = this->mesh_ptr->cells.nb();

        for (index_t c_i = 0; c_i < nb_cell; ++c_i)
        {
            std::array<double, 3> min_v;
            min_v.fill(std::numeric_limits<float>::infinity());

            std::array<double, 3> max_v;
            max_v.fill(-std::numeric_limits<float>::infinity());
            index_t nb_v = this->mesh_ptr->cells.nb_vertices(c_i);
            for (index_t v_i = 0;
                 v_i < nb_v; ++v_i)
            {
                index_t g_index = this->mesh_ptr->cells.vertex(c_i, v_i);
                auto v_pos = this->mesh_ptr->vertices.point(g_index);

                for (index_t i = 0; i < 3; i++)
                {
                    if (v_pos[i] < min_v[i])
                    {
                        min_v[i] = v_pos[i];
                    }

                    if (v_pos[i] > max_v[i])
                    {
                        max_v[i] = v_pos[i];
                    }
                }
            }

            c_range bx(pt(min_v[0], min_v[1], min_v[2]),
                       pt(max_v[0], max_v[1], max_v[2]));

            this->mesh_tree.insert(std::make_pair(bx, c_i));
        }
    }

    vector<value> MeshHeightSlicer::get_layer(double z)
    {
        vector<value> results;
        auto bds = this->mesh_tree.bounds();
        auto min_v = bds.min_corner();
        auto max_v = bds.max_corner();

        c_range bbx(pt(min_v.get<0>(), min_v.get<1>(), z),
                    pt(max_v.get<0>(), max_v.get<1>(), z));

        this->mesh_tree.query(bgi::intersects(bbx), std::back_inserter(results));
        Mesh m_out;
        clip_cell(results, m_out, z);

        return results;
    }

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
}