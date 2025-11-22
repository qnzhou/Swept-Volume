//
//  post_processing.h
//  general_sweep
//
//  Created by Yiwen Ju on 5/13/25.
//

#ifndef post_processing_h
#define post_processing_h
#include <ankerl/unordered_dense.h>
#include <arrangement/Arrangement.h>
#include <lagrange/SurfaceMesh.h>
#include <lagrange/combine_meshes.h>
#include <lagrange/compute_area.h>
#include <lagrange/separate_by_components.h>
#include <lagrange/views.h>

#include <algorithm>

template <typename Scalar, typename Index>
lagrange::SurfaceMesh<Scalar, Index> compute_swept_volume_from_envelope(
    const lagrange::SurfaceMesh<Scalar, Index>& envelope) {
    using Point = Eigen::Matrix<Scalar, 3, 1>;

    // Compute arrangement
    auto V = vertex_view(envelope).template cast<double>();
    auto F = facet_view(envelope).template cast<int>();
    auto T = attribute_vector_view<Scalar>(envelope, "time");
    arrangement::VectorI face_labels =
        Eigen::VectorXi::LinSpaced(F.rows(), 0, F.rows() - 1);
    auto engine =
        arrangement::Arrangement::create_mesh_arrangement(V, F, face_labels);
    engine->run();

    lagrange::SurfaceMesh<Scalar, Index> output_mesh;

    const auto& cell_data = engine->get_cells();  // (#facets x 2) array
    const auto& patches = engine->get_patches();  // list of facet indices
    const auto& parent_facets =
        engine->get_out_face_labels();  // size = #facets
    const auto& winding_number =
        engine->get_winding_number();  // (#facets x 2) array

    const auto& out_vertices = engine->get_vertices();
    const auto& arrangement_faces = engine->get_faces();
    size_t num_cells = engine->get_num_cells();
    size_t num_patches = engine->get_num_patches();
    size_t num_facets = arrangement_faces.rows();
    assert(patches.size() == num_facets);

    // Build cell adjacency graph and compute cell volumes
    constexpr Scalar vol_threshold = 1e-5;
    std::vector<ankerl::unordered_dense::set<Index>> cell_graph(num_cells);
    std::vector<Scalar> cell_volumes(num_cells, 0);
    for (size_t fid = 0; fid < num_facets; fid++) {
        Index c0 = static_cast<Index>(
            cell_data(patches[fid], 0));  // Cell on the positive side
        Index c1 = static_cast<Index>(
            cell_data(patches[fid], 1));  // Cell on the negative side

        cell_graph[c0].insert(c1);
        cell_graph[c1].insert(c0);

        Eigen::Matrix<int, 1, 3> f = arrangement_faces.row(fid);
        Point p0 = out_vertices.row(f[0]).template cast<Scalar>();
        Point p1 = out_vertices.row(f[1]).template cast<Scalar>();
        Point p2 = out_vertices.row(f[2]).template cast<Scalar>();
        Scalar vol = p0.dot((p1).cross(p2));

        cell_volumes[c0] += vol / 6;
        cell_volumes[c1] -= vol / 6;
    }
    std::vector<bool> cell_is_isolated(num_cells, false);
    for (size_t cid = 0; cid < num_cells; cid++) {
        if (cell_graph[cid].size() <= 1) {
            // Cell is completely embedded in another cell
            // It can be safely removed without affecting other cells and their
            // winding numbers
            cell_is_isolated[cid] = true;
        }
    }

    output_mesh.add_vertices(
        out_vertices.rows(),
        {out_vertices.data(), static_cast<size_t>(out_vertices.size())});

    size_t num_valid_facets = 0;
    std::vector<int8_t> valid(num_facets, 0);
    for (size_t fid = 0; fid < num_facets; fid++) {
        int w0 = winding_number(fid, 0);  // Winding number on the positive side
        int w1 = winding_number(fid, 1);  // Winding number on the negative side
        int c0 = cell_data(patches[fid], 0);  // Cell on the positive side
        int c1 = cell_data(patches[fid], 1);  // Cell on the negative side

        if (cell_is_isolated[c0] && std::abs(cell_volumes[c0]) < vol_threshold)
            continue;
        if (cell_is_isolated[c1] && std::abs(cell_volumes[c1]) < vol_threshold)
            continue;

        if (w0 == 0 && w1 != 0) {
            valid[fid] = 1;
            num_valid_facets++;
        } else if (w1 == 0 && w0 != 0) {
            valid[fid] = -1;
            num_valid_facets++;
        }
    }

    output_mesh.add_triangles(num_valid_facets);

    arrangement::MatrixIr out_facets(num_valid_facets, 3);
    Index count = 0;
    for (size_t fid = 0; fid < num_facets; fid++) {
        if (valid[fid] == 0) continue;
        auto f = output_mesh.ref_facet_vertices(count++);
        if (valid[fid] == 1) {
            f[0] = arrangement_faces(fid, 0);
            f[1] = arrangement_faces(fid, 1);
            f[2] = arrangement_faces(fid, 2);
        } else {
            f[0] = arrangement_faces(fid, 2);
            f[1] = arrangement_faces(fid, 1);
            f[2] = arrangement_faces(fid, 0);
        }
    }

    auto interpolate_time = [&](Index fid, Point p) {
        Eigen::Matrix<int, 1, 3> f = F.row(fid);
        Point p0 = V.row(f[0]);
        Point p1 = V.row(f[1]);
        Point p2 = V.row(f[2]);
        auto t0 = T[f[0]];
        auto t1 = T[f[1]];
        auto t2 = T[f[2]];

        Scalar b0 = ((p1 - p).cross(p2 - p)).norm();
        Scalar b1 = ((p2 - p).cross(p0 - p)).norm();
        Scalar b2 = ((p0 - p).cross(p1 - p)).norm();
        return (t0 * b0 + t1 * b1 + t2 * b2) / (b0 + b1 + b2);
    };

    output_mesh.template create_attribute<Scalar>(
        "time", lagrange::AttributeElement::Corner,
        lagrange::AttributeUsage::Scalar, 1);
    auto time_values = attribute_vector_ref<Scalar>(output_mesh, "time");
    count = 0;
    // TODO: parallalize this loop
    for (size_t fid = 0; fid < num_facets; fid++) {
        if (valid[fid] == 0) continue;
        auto idx = count++;
        auto parent_fid = parent_facets(fid);
        Index v0 = static_cast<Index>(arrangement_faces(fid, 0));
        Index v1 = static_cast<Index>(arrangement_faces(fid, 1));
        Index v2 = static_cast<Index>(arrangement_faces(fid, 2));
        Point p0(out_vertices.row(v0).template cast<Scalar>());
        Point p1(out_vertices.row(v1).template cast<Scalar>());
        Point p2(out_vertices.row(v2).template cast<Scalar>());
        Scalar t0 = interpolate_time(parent_fid, p0);
        Scalar t1 = interpolate_time(parent_fid, p1);
        Scalar t2 = interpolate_time(parent_fid, p2);
        if (valid[fid] == 1) {
            time_values[3 * idx + 0] = t0;
            time_values[3 * idx + 1] = t1;
            time_values[3 * idx + 2] = t2;
        } else {
            time_values[3 * idx + 0] = t2;
            time_values[3 * idx + 1] = t1;
            time_values[3 * idx + 2] = t0;
        }
    }

    return output_mesh;
}

#endif /* post_processing_h */
