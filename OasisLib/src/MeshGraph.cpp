#include "MeshGraph.h"
#include <geogram/basic/common.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/basic/logger.h>

#include <boost/geometry/algorithms/convex_hull.hpp>

namespace OasisLib
{

    MeshHeightSlicer::MeshHeightSlicer(sh_mesh_ptr m_ptr)
    {
        Logger::div("Init MeshHeightSlicer");
        this->mesh_ptr = m_ptr;

        Logger::out("MeshHeightSlicer") << "Input mesh has " \
                << this->mesh_ptr-> cells.nb() << \
                " cells," << \
                this->mesh_ptr->vertices.nb() <<
                " vertices" <<std::endl;

        this->initialize_rtree();
        this->nb_nodes = this->mesh_tree.size();
    }

    void MeshHeightSlicer::initialize_rtree()
    {
        index_t nb_c = this->mesh_ptr->cells.nb();
        for (index_t c = 0; c < nb_c; ++c)
        {
            index_t nb_e = this->mesh_ptr->cells.nb_edges(c);

            for (auto e_i = 0; e_i < nb_e; ++e_i)
            {
                pt min_p;
                pt max_p;
                auto p1_id = this->mesh_ptr->cells.edge_vertex(c, e_i, 0);
                auto p2_id = this->mesh_ptr->cells.edge_vertex(c, e_i, 1);

                auto p_1_pos = this->mesh_ptr->vertices.point(p1_id);
                auto p_2_pos = this->mesh_ptr->vertices.point(p2_id);

                for (int i = 0; i < 3; ++i)
                {
                    if (p_1_pos[i] < p_2_pos[i])
                    {
                        min_p[i] = p_1_pos[i];
                        max_p[i] = p_2_pos[i];
                    }
                    else
                    {
                        min_p[i] = p_2_pos[i];
                        max_p[i] = p_1_pos[i];
                    }
                }
                c_range bx(min_p, max_p);
                auto id = std::make_pair(c, e_i);
                this->mesh_tree.insert(std::make_pair(bx, id));
            }
        }
    }

    id_map MeshHeightSlicer::get_layer(Mesh &m_out, double z)
    {
        std::vector<value> results;
        auto bds = this->mesh_tree.bounds();
        auto min_v = bds.min_corner();
        auto max_v = bds.max_corner();

        pt l_pt({min_v.get<0>(), min_v.get<1>(), z});
        pt m_pt({max_v.get<0>(), max_v.get<1>(), z});

        c_range bbx(l_pt, m_pt);

        this->mesh_tree.query(bgi::intersects(bbx), std::back_inserter(results));

        Logger::out("MeshHeightSlicer") << "computed " << \
                results.size() << " intersection" << std::endl;

        auto c2f = this->clip_cell(results, m_out, z);

        return c2f;
    }

    id_map MeshHeightSlicer::clip_cell(std::vector<value> &target_facet, Mesh &m_out, double z)
    {
        id_map cell2facet;
        std::map<index_t, mpt2d> pt_accumulator;
        for (auto it = target_facet.cbegin(); it != target_facet.cend(); ++it)
        {
            auto c_id = it->second.first;
            auto e_id = it->second.second;

            auto c_it = pt_accumulator.find(c_id);

            if (c_it == pt_accumulator.cend())
            {
                pt_accumulator.insert({c_id, mpt2d()});
            }

            auto g_p1 = this->mesh_ptr->cells.edge_vertex(c_id, e_id, 0);
            auto g_p2 = this->mesh_ptr->cells.edge_vertex(c_id, e_id, 1);

            auto p_1 = this->mesh_ptr->vertices.point(g_p1);
            auto p_2 = this->mesh_ptr->vertices.point(g_p2);

            auto c_1_sign = geo_sgn(p_1[2] - z);
            auto c_2_sign = geo_sgn(p_2[2] - z);

            if (c_1_sign != c_2_sign && (c_1_sign != 0 && c_2_sign != 0))
            {
                pt2d inter_p2d;
                double d = (z - p_1[2]) / (p_2[2] - p_1[2]);
                auto inter_p = p_1 + d * (p_2 - p_1);

                inter_p2d[0] = inter_p[0];
                inter_p2d[1] = inter_p[1];
                bg::append(pt_accumulator[c_id], inter_p2d);
            }
            else
            {
                if (c_1_sign == 0)
                {
                    pt2d inter_p2d;
                    inter_p2d[0] = p_1[0];
                    inter_p2d[1] = p_1[1];
                    bg::append(pt_accumulator[c_id], inter_p2d);
                }

                if (c_2_sign == 0)
                {
                    pt2d inter_p2d;
                    inter_p2d[0] = p_2[0];
                    inter_p2d[1] = p_2[1];
                    bg::append(pt_accumulator[c_id], inter_p2d);
                }
            }
        }

        for (auto c_it = pt_accumulator.cbegin(); c_it != pt_accumulator.cend(); ++c_it)
        {
            auto c_id = c_it->first;
            auto m_pt = c_it->second;

            contour2d hull;
            bg::convex_hull(m_pt, hull);

            std::vector<index_t> added_id;

            for (auto p_it = std::cbegin(hull); p_it != std::cend(hull); ++p_it)
            {
                pt a_v;
                auto hull_pt = *p_it;
                a_v[0] = hull_pt[0];
                a_v[1] = hull_pt[1];
                a_v[2] = z;
                auto a_id = m_out.vertices.create_vertex(a_v.data());
                added_id.push_back(a_id);
            }

            auto f_id = m_out.facets.create_polygon(added_id.size(),&added_id[0]);
            cell2facet.insert({f_id, c_id});
        }
        auto repair_mode = static_cast<MeshRepairMode>
                        (static_cast<int>(MeshRepairMode::MESH_REPAIR_COLOCATE) | 
                        static_cast<int>(MeshRepairMode::MESH_REPAIR_DUP_F)) ;
        mesh_repair(m_out, repair_mode);
        return cell2facet;
    }

}