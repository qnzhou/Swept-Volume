#include "adaptive_column_grid.h"
#include "col_gridgen.h"
#include "io.h"
#include "post_processing.h"

#include <generalized_sweep/generalized_sweep.h>
#include <mtet/grid.h>
#include <mtetcol/contour.h>
#include <mtetcol/simplicial_column.h>

#include <array>
#include <iostream>
#include <vector>

namespace sweep {

std::tuple<std::vector<Scalar>, std::vector<Index>,
           std::vector<std::vector<Scalar>>, std::vector<std::vector<Scalar>>>
refine_grid(const SpaceTimeFunction& f, mtet::MTetMesh& grid,
            const SweepOptions& options) {
    std::cout << "Start to generate the background grid..." << std::endl;

    vertExtrude vertexMap;
    insidenessMap insideMap;

    // TODO: what is these timers?
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
    SweepResult result;
    auto grid = mtet::generate_tet_grid(grid_spec.resolution,
                                        grid_spec.bbox_min, grid_spec.bbox_max);
    // TODO: investigate why saving and loading is necessary here???
    mtet::save_mesh("init.msh", grid);
    grid = mtet::load_mesh("init.msh");

    auto [verts, simps, time, values] = refine_grid(f, grid, options);
    std::function<std::span<double>(size_t)> time_func =
        [&](size_t index) -> std::span<double> { return time[index]; };
    std::function<std::span<double>(size_t)> values_func =
        [&](size_t index) -> std::span<double> { return values[index]; };
    mtetcol::SimplicialColumn<4> columns;
    columns.set_vertices(verts);
    columns.set_simplices(simps);
    columns.set_time_samples(time_func, values_func);

    constexpr Scalar iso_value = 0.0;
    auto contour = columns.extract_contour(iso_value, options.cyclic);
    if (!contour.is_manifold()) {
        throw std::runtime_error("ERROR: extracted contour is not manifold");
    }

    result.envelope = compute_envelope(f, contour, options);

    result.arrangement = compute_envelope_arrangement(result.envelope);

    result.sweep_surface =
        extract_sweep_surface_from_arrangement(result.arrangement);

    return result;
}

}  // namespace sweep
