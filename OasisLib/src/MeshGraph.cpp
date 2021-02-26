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
            index_t nb_v = this->mesh_ptr->cells.nb_vertices(c);
            pt min_p;
            pt max_p;

            for (int i = 0; i < 3; ++i)
            {
                min_p[i] = std::numeric_limits<double>::infinity();
                max_p[i] = -std::numeric_limits<double>::infinity();
            }

            for (auto v_i = 0; v_i < nb_v; ++v_i)
            {
                index_t g_v = this->mesh_ptr->cells.vertex(c, v_i);
                auto v_pos = this->mesh_ptr->vertices.point(g_v);

                for (int i = 0; i < 3; ++i)
                {
                    if (min_p[i] > v_pos[i])
                    {
                        min_p[i] = v_pos[i];
                    }

                    if (max_p[i] < v_pos[i])
                    {
                        max_p[i] = v_pos[i];
                    }
                }
            }

            c_range bx(min_p, max_p);
            this->mesh_tree.insert(std::make_pair(bx, c));
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

        for (auto it = target_facet.cbegin(); it != target_facet.cend(); ++it)
        {
            auto cell_id = it->second;
            vector<index_t> poly;

            index_t inter_f;
            auto nb_f_c = this->mesh_ptr->cells.nb_facets(cell_id);
            for (auto f = 0; f < nb_f_c; ++f)
            {
                inter_f = f;
                auto intersect_f = IsFacet(cell_id, f, z);
                if (intersect_f)
                {
                    break;
                }
            }
            auto n_f = this->mesh_ptr->cells.facet(cell_id, inter_f);
            auto n_f_start = n_f;
            do {
                auto c_begin = this->mesh_ptr->facets.corners_begin(n_f);
                auto c_end = this->mesh_ptr->facets.corners_end(n_f);

                int count = 0;
                auto inter_c = 0;
                for (auto c_i = c_begin; c_i < c_end; ++c_i)
                {
                    pt inter_p;
                    auto c_i_1 = this->mesh_ptr->facets.next_corner_around_facet(n_f, c_i);
                    auto intersection = IsIntersection(c_i, c_i_1, inter_p, z);

                    if (intersection == 2)
                    {
                        auto g_p1 = this->mesh_ptr->facet_corners.vertex(c_i);
                        auto g_p2 = this->mesh_ptr->facet_corners.vertex(c_i_1);

                        auto p_1 = this->mesh_ptr->vertices.point(g_p1);
                        auto p_2 = this->mesh_ptr->vertices.point(g_p2);
                        auto v_id_1 = m_out.vertices.create_vertex(p_1.data());
                        auto v_id_2 = m_out.vertices.create_vertex(p_2.data());
                        poly.push_back(v_id_1);
                        poly.push_back(v_id_2);
                        count = 2;
                    }

                    if (intersection == 1)
                    {
                        auto v_id = m_out.vertices.create_vertex(inter_p.data());
                        poly.push_back(v_id);
                        count++;
                    }

                    if (count > 2)
                    {
                        inter_c = c_i;
                        break;
                    }  
                }

                n_f = this->mesh_ptr->facet_corners.adjacent_facet(inter_c);
            }while(n_f != n_f_start);

            m_out.facets.create_polygon(poly);
        }
    }

    int MeshHeightSlicer::IsIntersection(index_t c_1, index_t c_2, pt &inter_p, double z)
    {
        int intersection = 0;
        auto g_p1 = this->mesh_ptr->facet_corners.vertex(c_1);
        auto g_p2 = this->mesh_ptr->facet_corners.vertex(c_2);

        auto p_1 = this->mesh_ptr->vertices.point(g_p1);
        auto p_2 = this->mesh_ptr->vertices.point(g_p2);

        auto c_1_sign = geo_sgn(p_1[2] - z);
        auto c_2_sign = geo_sgn(p_2[2] - z);

        if (c_1_sign != c_2_sign && (c_1_sign != 0 && c_2_sign != 0))
        {
            double d = (z - p_1[2]) / (p_2[2] - p_1[2]);
            inter_p = p_1 + d * (p_2 - p_1);
            intersection = 1;
            return intersection;
        }

        if (c_1_sign == 0 && c_2_sign == 0)
        {
            intersection = 2;
            return intersection;
        }

        if (c_1_sign == 0)
        {
            inter_p = p_1;
            intersection = 1;
            return intersection;
        }

        if (c_2_sign == 0)
        {
            inter_p = p_2;
            intersection = 1;
            return intersection;
        }

        return intersection;
    }

    bool MeshHeightSlicer::IsFacet(index_t c, index_t f, double z)
    {
        auto nb_v = this->mesh_ptr->cells.facet_nb_vertices(c, f);

        double min_p = std::numeric_limits<double>::infinity();
        double max_p = -std::numeric_limits<double>::infinity();

        for (auto v = 0; v < nb_v; ++v)
        {
            auto v_id = this->mesh_ptr->cells.facet_vertex(c, f, v);
            auto v_pos = this->mesh_ptr->vertices.point(v_id);

            if (min_p > v_pos[2])
            {
                min_p = v_pos[2];
            }

            if (max_p < v_pos[2])
            {
                max_p = v_pos[2];
            }
        }

        if (min_p <= z && max_p >= z)
        {
            return true;
        }

        return false;
    }

}