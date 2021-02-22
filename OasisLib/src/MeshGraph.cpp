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
        index_t nb_cell = this->mesh_ptr->facets.nb();

        for (index_t f = 0; f < this->mesh_ptr->facets.nb(); ++f)
        {
            std::array<double, 3> min_v;
            min_v.fill(std::numeric_limits<float>::infinity());

            std::array<double, 3> max_v;
            max_v.fill(-std::numeric_limits<float>::infinity());

            for (index_t c = this->mesh_ptr->facets.corners_begin(f);
                 c < this->mesh_ptr->facets.corners_end(f); ++c)
            {
                index_t g_index = this->mesh_ptr->facet_corners.vertex(c);
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

            this->mesh_tree.insert(std::make_pair(bx, f));
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
        return results;
    }

    void MeshHeightSlicer::clip_cell(vector<value> &target_facet, Mesh &m_out, cell_id_map &id_map, double z)
    {
        // bisector
        for (auto it = target_facet.cbegin(); it != target_facet.cend(); ++it)
        {
            auto pair = *it;
            auto f = pair.second;
            //Polygon result;
            //clip_poly_by_plane(poly, result, this->mesh_ptr, p0.data(), p1.data());
        }
    }

    void clip_poly_by_plane(index_t f, const Mesh *mesh, double z)
    {

        vecng<3, double> p0{0.0, 0.0, z};
        vecng<3, double> n{0.0, 0.0, 1.0};

        auto c_begin = mesh ->facets.corners_begin(f);
        auto c_end = mesh ->facets.corners_end(f);
        std::pair<vecng<double, 3>,vecng<double, 3>> edge; 
        int c = 0;
        for(auto k = c_begin; k < c_end; ++k)
        {
            if (c > 1)
            {
                break;
            }

            index_t current_v_idx = mesh->facet_corners.vertex(k);
            auto c_pos = mesh->vertices.point(current_v_idx);

            index_t next_corner = mesh->facets.next_corner_around_facet(f, k);
            index_t next_v_idx = mesh->facet_corners.vertex(next_corner);
            auto n_pos = mesh->vertices.point(next_v_idx);

            auto is_intersect = PCK::dot_compare_3d(n.data(), n_pos.data(), c_pos.data());
            auto is_parrallel = PCK::dot_compare_3d(n.data(), p0.data(), v_pos.data());

            // if sign not identical then intersection
            if (is_intersect != 0)
            {
                auto inter_p = c_pos + dot((p0 - c_pos), n)/dot((n_pos - c_pos),n)*(n_pos - c_pos);
                if(c == 0)
                {
                    edge.first = inter_p;
                    c++;
                }else
                {
                    edge.second = inter_p;
                    c++
                }

            }else{
                if (is_parrallel == 0 && c == 0)
                {
                    edge.first = c_pos;
                    edge.second = n_pos;
                    c = 2
                }
            }
        }
    }
}