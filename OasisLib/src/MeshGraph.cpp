#include "MeshGraph.h"
#include <geogram/basic/common.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/attributes.h>

#include <boost/geometry/algorithms/convex_hull.hpp>

namespace OasisLib
{

    MeshHeightSlicer::MeshHeightSlicer(sh_mesh_ptr m_ptr)
    {
        Logger::div("Init MeshHeightSlicer");
        this->mesh_ptr = m_ptr;
        Logger::out("MeshHeightSlicer") << "Input mesh has " \
                                        << this->mesh_ptr->facets.nb() << " facets," << \
                                        this->mesh_ptr->vertices.nb() << " vertices" \
                                        << std::endl;

        check_mesh();
        this->initialize_rtree();
        this->nb_nodes = this->mesh_tree.size();
    }

    void MeshHeightSlicer::check_mesh()
    {
        try
        {
            if (this->mesh_ptr->cells.nb() == 0)
            {
                Attribute<int> cell_id(this->mesh_ptr->facets.attributes(), "cell_id");
                auto c_vector = cell_id.get_vector();
                std::set<int> counter(c_vector.cbegin(), c_vector.cend());
                Logger::out("MeshHeightSlicer") << counter.size() << " cells" << std::endl;

                if (counter.size() == 0)
                {
                    throw 0;
                }
            }
        }
        catch (int e)
        {
            Logger::err("MeshHeightSlicer") << "input mesh has no cells" << std::endl;
        }
    }

    void MeshHeightSlicer::initialize_rtree()
    {
        Attribute<int> cell_id(this->mesh_ptr->facets.attributes(), "cell_id");
        index_t nb_f = this->mesh_ptr->facets.nb();
        for (index_t f = 0; f < nb_f; ++f)
        {
            auto c = cell_id[f];
            auto c_begin = this -> mesh_ptr -> facets.corners_begin(f);
            auto c_end = this -> mesh_ptr -> facets.corners_end(f);
            for (auto c_i = c_begin; c_i != c_end; ++c_i)
            {
                pt min_p;
                pt max_p;
                auto p1_id = this->mesh_ptr->facet_corners.vertex(c_i);
                auto c_i_1 = this->mesh_ptr->facets.next_corner_around_facet(f, c_i);
                auto p2_id = this->mesh_ptr->facet_corners.vertex(c_i_1);

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
                std::array<index_t, 3> id = {(index_t)c, p1_id, p2_id};
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
        double rel_z = min_v.get<2>() + z;
        pt l_pt({min_v.get<0>(), min_v.get<1>(), rel_z - 0.1});
        pt m_pt({max_v.get<0>(), max_v.get<1>(), rel_z + 0.1});

        c_range bbx(l_pt, m_pt);

        this->mesh_tree.query(bgi::intersects(bbx), std::back_inserter(results));

        Logger::out("MeshHeightSlicer") << "computed " << results.size() << " intersection" << std::endl;

        auto c2f = this->clip_cell(results, m_out, rel_z);

        return c2f;
    }

    id_map MeshHeightSlicer::clip_cell(std::vector<value> &target_facet, Mesh &m_out, double z)
    {
        id_map cell2facet;
        std::map<index_t, mpt2d> pt_accumulator;
        for (auto it = target_facet.cbegin(); it != target_facet.cend(); ++it)
        {
            auto box = it -> first;
            auto id = it->second;
            auto c_id = id[0];

            auto c_it = pt_accumulator.find(c_id);

            if (c_it == pt_accumulator.cend())
            {
                pt_accumulator.insert({c_id, mpt2d()});
            }

            auto p_1 = this->mesh_ptr->vertices.point(id[1]);
            auto p_2 = this->mesh_ptr->vertices.point(id[2]);

            auto p_1_diff = p_1[2] - z;
            auto p_2_diff = p_2[2] - z;

            auto c_1_sign = geo_sgn(p_1_diff);
            auto c_2_sign = geo_sgn(p_2_diff);

            if (c_1_sign != c_2_sign && (abs(p_1_diff) >= 1e-4 && abs(p_1_diff) >= 1e-4))
            {
                pt2d inter_p2d;
                double d = (z - p_1[2]) / (p_2[2]- p_1[2]);
                auto inter_p = p_1 + d * (p_2 - p_1);

                inter_p2d[0] = inter_p[0];
                inter_p2d[1] = inter_p[1];
                bg::append(pt_accumulator[c_id], inter_p2d);
            }
            else
            {
                if (abs(p_1_diff) < 1e-4)
                {
                    pt2d inter_p2d;
                    inter_p2d[0] = p_1[0];
                    inter_p2d[1] = p_1[1];
                    bg::append(pt_accumulator[c_id], inter_p2d);
                }

                if (abs(p_1_diff) < 1e-4)
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

            auto count = bg::num_points(m_pt);
            if(count >= 3){
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

                auto f_id = m_out.facets.create_polygon(added_id.size(), &added_id[0]);
                cell2facet.insert({f_id, c_id});
            }
        }
        Logger::out("MeshHeightSlicer") << "outputs contains : " << cell2facet.size() << " facet" << std::endl;
        mesh_repair(m_out, MESH_REPAIR_COLOCATE, 1e-4);
        m_out.facets.connect();
        return cell2facet;
    }

}