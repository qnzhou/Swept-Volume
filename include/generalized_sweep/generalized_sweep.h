#pragma once

#include <lagrange/SurfaceMesh.h>

#include <Eigen/Core>
#include <array>
#include <functional>
#include <limits>

namespace generalized_sweep {

using Scalar = double;
using Index = uint32_t;

///
/// Type alias for a space-time implicit function.
///
/// The function takes a 4D point (x, y, z, t) as input and returns a pair:
///     - The first element is the scalar value of the implicit function at that
///     point.
///     - The second element is a 4D gradient vector (∂f/∂x, ∂f/∂y, ∂f/∂z,
///     ∂f/∂t) at that point.
///
using SpaceTimeFunction =
    std::function<std::pair<double, Eigen::RowVector4d>(Eigen::RowVector4d)>;

struct GridSpec {
    std::array<size_t, 3> resolution = {4, 4, 4};
    std::array<float, 3> bbox_min = {-0.2f, -0.2f, -0.2f};
    std::array<float, 3> bbox_max = {1.2f, 1.2f, 1.2f};
};

struct SweepResult {
    lagrange::SurfaceMesh<Scalar, Index> envelope;
    lagrange::SurfaceMesh<Scalar, Index> arrangement;
    lagrange::SurfaceMesh<Scalar, Index> sweep_surface;
};

struct SweepOptions {
    Scalar isovalue = 0.0;
    Scalar epsilon_env = 5e-4;
    Scalar epsilon_sil = 5e-3;
    int max_split = std::numeric_limits<int>::max();
    bool with_insideness_check = true;
    bool with_snapping = true;
    bool cyclic = false;
};

SweepResult generalized_sweep(const SpaceTimeFunction& f,
                              GridSpec grid_spec = {},
                              SweepOptions options = {});

};  // namespace generalized_sweep
