#pragma once

#include <lagrange/SurfaceMesh.h>

#include <Eigen/Core>
#include <array>
#include <functional>
#include <limits>

namespace sweep {

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

///
/// Specification for the 3D spatial grid used in the sweep.
///
struct GridSpec {
    /// Resolution of the grid in the x, y, z directions.
    std::array<size_t, 3> resolution = {4, 4, 4};

    /// Axis-aligned bounding box minimum corner.
    std::array<float, 3> bbox_min = {-0.2f, -0.2f, -0.2f};

    /// Axis-aligned bounding box maximum corner.
    std::array<float, 3> bbox_max = {1.2f, 1.2f, 1.2f};
};

///
/// Result of the generalized sweep operation.
///
struct SweepResult {
    /// The (lifted) envelope mesh.
    ///
    /// The spatial coordinates are stored as vertices of the mesh.
    /// The temporal coordinate is stored as a per-vertex attribute named
    /// "time".
    lagrange::SurfaceMesh<Scalar, Index> envelope;

    /// The spatial arrangement of the envelope.
    ///
    /// The temporal coordinate is stored as a per-corner attribute named
    /// "time".
    ///
    /// The per-facet attribute "valid" indicates the set of facets that are
    /// part of the final sweep surface:
    ///   - 0: the facet is not part of the sweep surface
    ///   - 1: the facet is part of the sweep surface with original orientation
    ///   - -1: the facet is part of the sweep surface with flipped orientation
    ///
    /// The per-edge attribute "is_feature" indicates which edge is part of
    /// feature curves on the sweep surface.
    lagrange::SurfaceMesh<Scalar, Index> arrangement;

    /// The final sweep surface mesh extracted from the arrangement.
    ///
    /// The per-corner attribute "time" stores the temporal coordinate.
    lagrange::SurfaceMesh<Scalar, Index> sweep_surface;
};

///
/// Options for controlling the generalized sweep operation.
///
struct SweepOptions {
    /// Tolerance for envelope grid refinement.
    Scalar epsilon_env = 5e-4;

    /// Tolerance for silhouette grid refinement.
    Scalar epsilon_sil = 5e-3;

    /// Maximum number of splits allowed during grid refinement.
    int max_split = std::numeric_limits<int>::max();

    /// Whether to perform insideness checks during grid refinement.
    /// If true, the algorithm will stop refinement early once it detects a cell
    /// is inside of the swept volume.
    bool with_insideness_check = true;

    /// Whether to enable vertex snapping during isocontouring.
    bool with_snapping = true;

    /// Whether the trajectory is cyclic.
    ///
    /// This feature is experimental and is not fully supported.
    bool cyclic = false;
};

///
/// Compute the generalized sweep surface.
///
/// @param f         The space-time implicit function defining the sweep.
/// @param grid_spec Specification of the initial 3D spatial grid.
/// @param options   Options controlling the sweep operation.
///
/// @return          The result of the generalized sweep, including the envelope
/// mesh,
///                  the arrangement mesh, and the final sweep surface.
///
SweepResult generalized_sweep(const SpaceTimeFunction& f,
                              GridSpec grid_spec = {},
                              SweepOptions options = {});

};  // namespace sweep
