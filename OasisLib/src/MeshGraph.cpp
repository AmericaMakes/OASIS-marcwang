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
        index_t nb_f = this->mesh_ptr->facets.nb();
        for (index_t f_i = 0; f_i < nb_f; ++f_i)
        {
            for (auto c_i = this->mesh_ptr->facets.corners_begin(f_i);
                 c_i < this->mesh_ptr->facets.corners_end(f_i); ++c_i)
            {
                MeshHalfedges::Halfedge hl(f_i, c_i);

                auto p_1_id = this->mesh_ptr->facet_corners.vertex(c_i);
                auto c_i_1 = this->mesh_ptr->facets.next_corner_around_facet(f_i, c_i);

                auto p_2_id = this->mesh_ptr->facet_corners.vertex(c_i_1);

                auto p_1 = this->mesh_ptr->vertices.point(p_1_id);
                auto p_2 = this->mesh_ptr->vertices.point(p_2_id);

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
                auto hl_id = std::make_pair(hl.facet, hl.corner);
                this->mesh_tree.insert(std::make_pair(bx, hl_id));
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
        Mesh m_out;
        clip_cell(results, m_out, z);

        return results;
    }

    void MeshHeightSlicer::clip_cell(vector<value> &target_facet, Mesh &m_out, double z)
    {
        MeshHalfedges mhl(*this->mesh_ptr);

        std::map<fc_p, pt> processed;
        for (auto it = target_facet.cbegin(); it != target_facet.cend(); ++it)
        {
            auto pair = *it;
            auto id_pair = pair.second;

            auto f_it = processed.find(id_pair);

            MeshHalfedges::Halfedge c_hl(id_pair.first, id_pair.second);

            if (f_it == processed.cend())
            {
                MeshHalfedges::Halfedge n_hl(c_hl.facet, c_hl.corner);
                mhl.move_to_next_around_facet(n_hl);

                auto p_1_id = this->mesh_ptr->facet_corners.vertex(c_hl.corner);
                auto p_1 = this->mesh_ptr->vertices.point(p_1_id);

                auto c_1_sign = geo_sgn(p_1[2] - z);
                int count = 0;

                while (count < 2)
                {

                    auto p_2_id = this->mesh_ptr->facet_corners.vertex(n_hl.corner);
                    auto p_2 = this->mesh_ptr->vertices.point(p_2_id);

                    auto c_2_sign = geo_sgn(p_2[2] - z);

                    MeshHalfedges::Halfedge opp_hl(c_hl.facet, c_hl.corner);
                    mhl.move_to_opposite(opp_hl);

                    if (c_1_sign != c_2_sign && c_1_sign != 0)
                    {
                        double d = (z - p_1[2]) / (p_2[2] - p_1[2]);
                        auto inter_p = p_1 + d * (p_2 - p_1);
                        processed.insert({{c_hl.facet, c_hl.corner}, inter_p});
                        processed.insert({{opp_hl.facet, opp_hl.corner}, inter_p});
                        count++;
                    }

                    if (c_1_sign == 0)
                    {
                        processed.insert({{c_hl.facet, c_hl.corner}, p_1});
                        processed.insert({{opp_hl.facet, opp_hl.corner}, p_1});
                        count++;
                    }

                    if (processed.find({n_hl.facet, n_hl.corner}) != processed.cend())
                    {
                        count++;
                    }
                    c_hl = n_hl;
                    mhl.move_to_next_around_facet(n_hl);
                    p_1 = p_2;
                    c_1_sign = c_2_sign;
                }
            }
        }
    }
}