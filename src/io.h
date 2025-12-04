//
//  io.h
//  adaptive_colume_grid
//
//  Created by Yiwen Ju on 12/3/24.
//

#ifndef io_h
#define io_h


#include "adaptive_column_grid.h"
#include <mshio/mshio.h>

#include <fstream>

#include <lagrange/SurfaceMesh.h>
#include <lagrange/views.h>

void convert_4d_grid_col(mtet::MTetMesh grid,
                         vertExtrude vertexMap,
                         std::vector<std::array<double, 3>> &verts,
                         std::vector<std::array<size_t, 4>> &simps,
                         std::vector<std::vector<double>> &time,
                         std::vector<std::vector<double>> &values);

void convert_4d_grid_mtetcol(mtet::MTetMesh grid,
                             vertExtrude vertexMap,
                             std::vector<double> &verts,
                             std::vector<uint32_t> &simps,
                             std::vector<std::vector<double>> &time,
                             std::vector<std::vector<double>> &values,
                             bool cyclic);

mshio::MshSpec generate_spec(
                             const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F,
                             TimeMap timeMap);

Eigen::VectorXd propagate_labels_bfs(
                                     const Eigen::MatrixXd& V,   // n x 3
                                     const Eigen::MatrixXi& F,   // m x 3
                                     const TimeMap& timeMap);

void backfill_timeMap_from_labels(
                                  const Eigen::MatrixXd& V,
                                  const Eigen::VectorXd& L,
                                  TimeMap& timeMap);

void save_grid_for_mathematica(
        std::string_view filename,
        mtet::MTetMesh grid,
        vertExtrude vertexMap);

template <typename Scalar, typename Index>
void save_features(std::string_view filename, lagrange::SurfaceMesh<Scalar, Index>& arrangement)
{
    std::ofstream fout(filename.data());

    assert(arrangement.has_attribute("is_feature"));
    auto vertex_view = lagrange::vertex_view(arrangement);
    Index num_vertices = arrangement.get_num_vertices();
    for (Index vid=0; vid < num_vertices; vid++) {
        fout << "v " << vertex_view(vid,0) << " "
             << vertex_view(vid,1) << " "
             << vertex_view(vid,2) << std::endl;
    }

    auto is_feature = lagrange::attribute_vector_view<int8_t>(arrangement, "is_feature");
    Index num_edges = arrangement.get_num_edges();
    for (Index eid=0; eid < num_edges; eid++) {
        if (is_feature(eid)) {
            auto [v0, v1] = arrangement.get_edge_vertices(eid);
            fout << "l "<< (v0 + 1) << " " << (v1 + 1) << std::endl;
        }
    }
}

#endif /* io_h */
