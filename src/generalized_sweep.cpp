#include <ankerl/unordered_dense.h>
#include <mtet/grid.h>
#include <mtetcol/contour.h>
#include <mtetcol/simplicial_column.h>
#include <nanothread/nanothread.h>
#include <sweep/generalized_sweep.h>
#include <sweep/logger.h>

#include <array>
#include <chrono>
#include <iostream>
#include <vector>

#include "adaptive_column_grid.h"
#include "col_gridgen.h"
#include "io.h"
#include "post_processing.h"

namespace sweep {

std::tuple<std::vector<Scalar>, std::vector<Index>,
           std::vector<std::vector<Scalar>>, std::vector<std::vector<Scalar>>>
refine_grid(const SpaceTimeFunction& f, mtet::MTetMesh& grid,
            const SweepOptions& options) {
    logger().info("Adaptively refine the background grid...");

    // TODO: investigate why saving and loading is necessary here???
    mtet::save_mesh("init.msh", grid);
    grid = mtet::load_mesh("init.msh");

    vertExtrude vertexMap;
    insidenessMap insideMap;

    // TODO: Clarify the purpose of these timers and whether they're still
    // needed.
    std::array<double, timer_amount> profileTimer{};
    std::array<size_t, timer_amount> profileCount{};
    spdlog::set_level(spdlog::level::off);
    if (!gridRefine(grid, vertexMap, insideMap, f, options.epsilon_env,
                    options.epsilon_sil, options.max_split,
                    options.with_insideness_check, profileTimer,
                    profileCount)) {
        throw std::runtime_error("ERROR: grid generation failed");
    };
    spdlog::set_level(spdlog::level::info);

    bool cyclic = options.cyclic;
    std::vector<mtetcol::Scalar> verts;
    std::vector<mtetcol::Index> simps;
    std::vector<std::vector<double>> time;
    std::vector<std::vector<double>> values;
    convert_4d_grid_mtetcol(grid, vertexMap, verts, simps, time, values,
                            cyclic);

    return {verts, simps, time, values};
}

std::tuple<std::vector<Scalar>, std::vector<Index>,
           std::vector<std::vector<Scalar>>, std::vector<std::vector<Scalar>>>
evaluate_grid(const SpaceTimeFunction& f, mtet::MTetMesh& grid,
              const SweepOptions& options) {
    size_t num_vertices = grid.get_num_vertices();
    size_t num_tets = grid.get_num_tets();

    std::vector<mtetcol::Scalar> verts(num_vertices * 3);
    std::vector<mtetcol::Index> simps(num_tets * 4);
    std::vector<std::vector<double>> time(num_vertices);
    std::vector<std::vector<double>> values(num_vertices);
    for (size_t i = 0; i < num_vertices; i++) {
        time[i].reserve(options.initial_time_samples);
        values[i].reserve(options.initial_time_samples);
    }

    // Extract vertices
    ankerl::unordered_dense::map<uint64_t, Index> vertex_map;
    vertex_map.reserve(num_vertices);
    grid.seq_foreach_vertex(
        [&](mtet::VertexId vid, std::span<const mtet::Scalar, 3> pos) {
            Index idx = static_cast<Index>(vertex_map.size());
            vertex_map[value_of(vid)] = idx;
            verts[idx * 3 + 0] = static_cast<mtetcol::Scalar>(pos[0]);
            verts[idx * 3 + 1] = static_cast<mtetcol::Scalar>(pos[1]);
            verts[idx * 3 + 2] = static_cast<mtetcol::Scalar>(pos[2]);
        });

    // Extract tets
    size_t tet_count = 0;
    grid.seq_foreach_tet([&](mtet::TetId tid,
                             std::span<const mtet::VertexId, 4> tet_verts) {
        for (size_t i = 0; i < 4; ++i) {
            simps[tet_count * 4 + i] =
                static_cast<mtetcol::Index>(vertex_map[value_of(tet_verts[i])]);
        }
        tet_count++;
    });

    // Extract time and time derivative
    namespace dr = drjit;
    dr::parallel_for(
        dr::blocked_range<size_t>(0, num_vertices, 1),
        [&](dr::blocked_range<size_t> idx_range) {
            for (size_t vid = idx_range.begin(); vid < idx_range.end(); vid++) {
                for (size_t tid = 0; tid < options.initial_time_samples;
                     tid++) {
                    double t = static_cast<double>(tid) /
                               (options.initial_time_samples - 1);
                    Eigen::RowVector4d eval_point;
                    eval_point << verts[vid * 3 + 0], verts[vid * 3 + 1],
                        verts[vid * 3 + 2], t;
                    auto eval = f(eval_point);
                    time[vid].push_back(t);
                    values[vid].push_back(
                        eval.second[3]);  // Value is time derivative.
                }
            }
        });

    return {verts, simps, time, values};
}

lagrange::SurfaceMesh<Scalar, Index> compute_envelope(
    const SpaceTimeFunction& f, mtetcol::Contour<4>& contour,
    const SweepOptions& options) {
    constexpr int dim = 4;
    size_t num_contour_vertices = contour.get_num_vertices();
    std::vector<double> function_values(num_contour_vertices);
    std::vector<double> gradient_values(num_contour_vertices * dim);
    for (size_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        auto pos_eval = f(Eigen::RowVector4d{pos[0], pos[1], pos[2], pos[3]});
        function_values[i] = pos_eval.first;
        gradient_values[dim * i] = pos_eval.second[0];
        gradient_values[dim * i + 1] = pos_eval.second[1];
        gradient_values[dim * i + 2] = pos_eval.second[2];
        gradient_values[dim * i + 3] = pos_eval.second[3];
    }

    // Extract isocontour
    auto isocontour = contour.isocontour(function_values, gradient_values,
                                         options.with_snapping);
    if (!isocontour.is_manifold()) {
        throw std::runtime_error("ERROR: extracted isocontour is not manifold");
    }
    if (isocontour.get_num_cycles() == 0) {
        throw std::runtime_error("ERROR: extracted isocontour has zero cycles");
    }
    isocontour.triangulate_cycles(true);
    lagrange::SurfaceMesh<Scalar, Index> envelope =
        isocontour_to_mesh<Scalar, Index>(isocontour);
    envelope.initialize_edges();

    return envelope;
}

SweepResult generalized_sweep(const SpaceTimeFunction& f, GridSpec grid_spec,
                              SweepOptions options) {
    auto init_grid_start = std::chrono::high_resolution_clock::now();
    SweepResult result;
    auto grid = mtet::generate_tet_grid(grid_spec.resolution,
                                        grid_spec.bbox_min, grid_spec.bbox_max);
    auto init_grid_end = std::chrono::high_resolution_clock::now();
    logger().info(
        "Initial grid generation time: {} seconds",
        std::chrono::duration<double>(init_grid_end - init_grid_start).count());

    std::vector<mtetcol::Scalar> verts;
    std::vector<mtetcol::Index> simps;
    std::vector<std::vector<double>> time;
    std::vector<std::vector<double>> values;

    if (options.with_adaptive_refinement) {
        auto refine_start = std::chrono::high_resolution_clock::now();
        std::tie(verts, simps, time, values) = refine_grid(f, grid, options);
        auto refine_end = std::chrono::high_resolution_clock::now();
        logger().info(
            "Grid refinement time: {} seconds",
            std::chrono::duration<double>(refine_end - refine_start).count());
    } else {
        auto evaluate_start = std::chrono::high_resolution_clock::now();
        std::tie(verts, simps, time, values) = evaluate_grid(f, grid, options);
        auto evaluate_end = std::chrono::high_resolution_clock::now();
        logger().info(
            "Grid evaluation time: {} seconds",
            std::chrono::duration<double>(evaluate_end - evaluate_start)
                .count());
    }

    std::function<std::span<double>(size_t)> time_func =
        [&](size_t index) -> std::span<double> { return time[index]; };
    std::function<std::span<double>(size_t)> values_func =
        [&](size_t index) -> std::span<double> { return values[index]; };
    mtetcol::SimplicialColumn<4> columns;
    columns.set_vertices(verts);
    columns.set_simplices(simps);
    columns.set_time_samples(time_func, values_func);

    constexpr Scalar iso_value = 0.0;

    auto silhouette_start = std::chrono::high_resolution_clock::now();
    auto contour = columns.extract_contour(iso_value, options.cyclic);
    auto silhouette_end = std::chrono::high_resolution_clock::now();
    logger().info(
        "Silhouette extraction time: {} seconds",
        std::chrono::duration<double>(silhouette_end - silhouette_start)
            .count());

    if (!contour.is_manifold()) {
        throw std::runtime_error("ERROR: extracted contour is not manifold");
    }

    auto envelope_start = std::chrono::high_resolution_clock::now();
    result.envelope = compute_envelope(f, contour, options);
    auto envelope_end = std::chrono::high_resolution_clock::now();
    logger().info(
        "Envelope computation time: {} seconds",
        std::chrono::duration<double>(envelope_end - envelope_start).count());

    auto arrangement_start = std::chrono::high_resolution_clock::now();
    result.arrangement =
        compute_envelope_arrangement(result.envelope, options.volume_threshold,
                                     options.face_count_threshold);
    auto arrangement_end = std::chrono::high_resolution_clock::now();
    logger().info(
        "Arrangement computation time: {} seconds",
        std::chrono::duration<double>(arrangement_end - arrangement_start)
            .count());

    auto sweep_surface_start = std::chrono::high_resolution_clock::now();
    result.sweep_surface =
        extract_sweep_surface_from_arrangement(result.arrangement);
    auto sweep_surface_end = std::chrono::high_resolution_clock::now();
    logger().info(
        "Sweep surface extraction time: {} seconds",
        std::chrono::duration<double>(sweep_surface_end - sweep_surface_start)
            .count());

    return result;
}

}  // namespace sweep
