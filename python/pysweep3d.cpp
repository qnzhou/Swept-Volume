#include <sweep/generalized_sweep.h>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

NB_MODULE(pysweep3d, m)
{
    using namespace nb::literals;
    using Scalar = sweep::Scalar;
    using Index = sweep::Index;

    m.doc() = "Python bindings for the generalized sweep library";

    // Bind GridSpec
    nb::class_<sweep::GridSpec>(
        m,
        "GridSpec",
        "Specification for the 3D spatial grid used in the sweep")
        .def(nb::init<>())
        .def_rw(
            "resolution",
            &sweep::GridSpec::resolution,
            "Resolution of the grid in the x, y, z directions")
        .def_rw("bbox_min", &sweep::GridSpec::bbox_min, "Axis-aligned bounding box minimum corner")
        .def_rw("bbox_max", &sweep::GridSpec::bbox_max, "Axis-aligned bounding box maximum corner");

    // Bind SweepOptions
    nb::class_<sweep::SweepOptions>(
        m,
        "SweepOptions",
        "Options for controlling the generalized sweep operation")
        .def(nb::init<>())
        .def_rw(
            "epsilon_env",
            &sweep::SweepOptions::epsilon_env,
            "Tolerance for envelope grid refinement")
        .def_rw(
            "epsilon_sil",
            &sweep::SweepOptions::epsilon_sil,
            "Tolerance for silhouette grid refinement")
        .def_rw(
            "max_split",
            &sweep::SweepOptions::max_split,
            "Maximum number of splits allowed during grid refinement")
        .def_rw(
            "with_insideness_check",
            &sweep::SweepOptions::with_insideness_check,
            "Whether to perform insideness checks during grid refinement")
        .def_rw(
            "with_snapping",
            &sweep::SweepOptions::with_snapping,
            "Whether to enable vertex snapping during isocontouring")
        .def_rw("cyclic", &sweep::SweepOptions::cyclic, "Whether the trajectory is cyclic")
        .def_rw(
            "volume_threshold",
            &sweep::SweepOptions::volume_threshold,
            "Minimum volume threshold for arrangement cell filtering")
        .def_rw(
            "face_count_threshold",
            &sweep::SweepOptions::face_count_threshold,
            "Minimum face count threshold for arrangement cell filtering")
        .def_rw(
            "with_adaptive_refinement",
            &sweep::SweepOptions::with_adaptive_refinement,
            "Adaptively refine the input grid based on the implicit function")
        .def_rw(
            "initial_time_samples",
            &sweep::SweepOptions::initial_time_samples,
            "Number of initial uniform time samples per spatial grid vertex")
        .def_rw(
            "min_tet_radius_ratio",
            &sweep::SweepOptions::min_tet_radius_ratio,
            "The minimum acceptable tetrahedron radius ratio during grid refinement")
        .def_rw(
            "min_tet_edge_length",
            &sweep::SweepOptions::min_tet_edge_length,
            "Minimum acceptable tetrahedron edge length during grid refinement");

    // Bind SweepResult
    nb::class_<sweep::SweepResult>(m, "SweepResult", "Result of the generalized sweep operation")
        .def("__init__", [](sweep::SweepResult& self) { new (&self) sweep::SweepResult(); })
        .def_rw("envelope", &sweep::SweepResult::envelope, "The (lifted) envelope mesh")
        .def_rw(
            "arrangement",
            &sweep::SweepResult::arrangement,
            "The spatial arrangement of the envelope")
        .def_rw(
            "sweep_surface",
            &sweep::SweepResult::sweep_surface,
            "The final sweep surface mesh extracted from the arrangement");

    // Bind the main generalized_sweep function
    m.def(
        "generalized_sweep",
        &sweep::generalized_sweep,
        "f"_a,
        "grid_spec"_a = sweep::GridSpec{},
        "options"_a = sweep::SweepOptions{},
        R"(
Compute the generalized sweep surface.

:param f: The space-time implicit function defining the sweep.
    Takes a 4D point (x, y, z, t) as input and returns a pair:
    - The scalar value of the implicit function at that point
    - A 4D gradient vector (∂f/∂x, ∂f/∂y, ∂f/∂z, ∂f/∂t) at that point
:type f: SpaceTimeFunction
:param grid_spec: Specification of the initial 3D spatial grid
:type grid_spec: GridSpec, optional
:param options: Options controlling the sweep operation
:type options: SweepOptions, optional
:return: The result of the generalized sweep, including the envelope mesh,
    the arrangement mesh, and the final sweep surface
:rtype: SweepResult
)");

    m.def(
        "generalized_sweep_from_config",
        &sweep::generalized_sweep_from_config,
        "function_file"_a,
        "config_file"_a,
        R"(
Compute the generalized sweep surface from configuration files.

:param function_file: Path to a file defining the space-time implicit function
:type function_file: str
:param config_file: Path to a configuration file specifying grid and options
:type config_file: str
:return: The result of the generalized sweep, including the envelope mesh,
    the arrangement mesh, and the final sweep surface
:rtype: SweepResult
)");
}
