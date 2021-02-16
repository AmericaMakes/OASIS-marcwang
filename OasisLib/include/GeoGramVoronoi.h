#pragma once
#ifndef GEOGRAMVORONOI_H
#define GEOGRAMVORONOI_H

int polyhedral_mesher(Mesh &M_in, Mesh &M_out,
                      int nb_points = 1000,
                      std::string simplify = "tets_voronoi_boundary",
                      double angle_threshold = 0.001,
                      index_t nb_iter_lloyd = 5,
                      index_t nb_iter_newton = 30,
                      bool tessallate_non_convex = false,
                      float poly_cell_shrinks = 0.0f,
                      bool generate_ids = true);

#endif