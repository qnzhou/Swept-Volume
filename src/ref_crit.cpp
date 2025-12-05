//
//  ref_crit.cpp
//  adaptive_column_grid
//
//  Created by Yiwen Ju on 12/10/24.
//
// #define no_inside
#include <ref_crit.h>

///Stores the 2d origin for the convex hull check happened in zero-crossing criteria.
std::array<double, 2> query_2d = {0.0, 0.0}; // X, Y

/// bezier ordinates of a 4D cubic simplex
const Eigen::Matrix<double, 35, 5, Eigen::RowMajor> ls{/// Vertices
                                                       {3, 0, 0, 0, 0},
                                                       {0, 3, 0, 0, 0},
                                                       {0, 0, 3, 0, 0},
                                                       {0, 0, 0, 3, 0},
                                                       {0, 0, 0, 0, 3},

                                                       /// Edges
                                                       {2, 1, 0, 0, 0},
                                                       {2, 0, 1, 0, 0},
                                                       {2, 0, 0, 1, 0},
                                                       {2, 0, 0, 0, 1},
                                                       {1, 2, 0, 0, 0},
                                                       {0, 2, 1, 0, 0},
                                                       {0, 2, 0, 1, 0},
                                                       {0, 2, 0, 0, 1},
                                                       {1, 0, 2, 0, 0},
                                                       {0, 1, 2, 0, 0},
                                                       {0, 0, 2, 1, 0},
                                                       {0, 0, 2, 0, 1},
                                                       {1, 0, 0, 2, 0},
                                                       {0, 1, 0, 2, 0},
                                                       {0, 0, 1, 2, 0},
                                                       {0, 0, 0, 2, 1},
                                                       {1, 0, 0, 0, 2},
                                                       {0, 1, 0, 0, 2},
                                                       {0, 0, 1, 0, 2},
                                                       {0, 0, 0, 1, 2},

                                                       /// Faces
                                                       {1, 1, 1, 0, 0},
                                                       {1, 1, 0, 1, 0},
                                                       {1, 1, 0, 0, 1},
                                                       {1, 0, 1, 1, 0},
                                                       {1, 0, 1, 0, 1},
                                                       {1, 0, 0, 1, 1},
                                                       {0, 1, 1, 1, 0},
                                                       {0, 1, 1, 0, 1},
                                                       {0, 1, 0, 1, 1},
                                                       {0, 0, 1, 1, 1}};

/// Hard coded parameters that take the derivative of a bezier simplex.
const std::array<std::array<size_t, 2>, 15> derMatrix{
    {{0, 8},
     {5, 27},
     {6, 29},
     {7, 30},
     {8, 21},
     {9, 12},
     {25, 32},
     {26, 33},
     {27, 22},
     {13, 16},
     {28, 34},
     {29, 23},
     {17, 20},
     {30, 24},
     {21, 4}}};

/// Hard coded parameters that elevate the simplex to cubic.
const std::array<std::array<size_t, 5>, 35> elevMatrix{
    {{0, 15, 15, 15, 15},  {15, 5, 15, 15, 15}, {15, 15, 9, 15, 15}, {15, 15, 15, 12, 15},
     {15, 15, 15, 15, 14}, {1, 0, 15, 15, 15},  {2, 15, 0, 15, 15},  {3, 15, 15, 0, 15},
     {4, 15, 15, 15, 0},   {5, 1, 15, 15, 15},  {15, 6, 5, 15, 15},  {15, 7, 15, 5, 15},
     {15, 8, 15, 15, 5},   {9, 15, 2, 15, 15},  {15, 9, 6, 15, 15},  {15, 15, 10, 9, 15},
     {15, 15, 11, 15, 9},  {12, 15, 15, 3, 15}, {15, 12, 15, 7, 15}, {15, 15, 12, 10, 15},
     {15, 15, 15, 13, 12}, {14, 15, 15, 15, 4}, {15, 14, 15, 15, 8}, {15, 15, 14, 15, 11},
     {15, 15, 15, 14, 13}, {6, 2, 1, 15, 15},   {7, 3, 15, 1, 15},   {8, 4, 15, 15, 1},
     {10, 15, 3, 2, 15},   {11, 15, 4, 15, 2},  {13, 15, 15, 4, 3},  {15, 10, 7, 6, 15},
     {15, 11, 8, 15, 6},   {15, 13, 15, 8, 7},  {15, 15, 13, 11, 10}}};

// Bézier stage-2 linear map  (outputs 35 control values from the 30-vector G)

/// returns a `bool` value that `true` represents positive and `false` represents negative of the input value `x`.
inline bool get_sign(const double x)
{
    return x > 0;
}

/// Hard coded indices of top and bottom tetrahedra of the 4D Bezier simplex control coordinates(35 indices in total).
const std::array<int, 16> topFIndices = {0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 20, 21, 23, 26};
const std::array<int, 16> botFIndices =
    {5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 26, 27, 28, 29};

/// Construct the 4D Bezier control points for a given 5-point (penta-linear) element in spacetime.
/// This function interpolates both position and gradient information across the 5 vertices
/// of a 4D simplex (5-cell) to generate 35 Bezier control points representing the 4D Bezier simplex.
/// The resulting Bezier simplex can then be used for zero crossing check and distance check.
///
/// @param[in] p1–p5        The 4D vertex positions (x, y, z, t) of the 5-cell. Each is an Eigen::RowVector4d.
/// @param[in] v1–v5        The scalar function values (e.g., implicit field or signed distance) associated with each vertex.
/// @param[in] g1–g5        The gradients of the scalar field at each vertex, stored as 4D row vectors corresponding to ∂f/∂x, ∂f/∂y, ∂f/∂z, ∂f/∂t.
/// @param[out] bezier      An array of 35 Bezier control points representing the 4D Bezier simplex of the 5-cell.
///                         These control points encode both positional and derivative information to ensure
///                         smooth interpolation across the simplex.
/// @param[out] inside      A conservative tag indicating whether all Bezier control points are negative,
///                         implying that the entire 5-cell is inside the swept volume (no sign change within the simplex).
///
/// @return                 A boolean value that indicates whether the Bezier simplex has a zero crossing using sweep function

[[gnu::always_inline]]
inline bool bezier4D(
    const Eigen::RowVector4d& p1,
    const Eigen::RowVector4d& p2,
    const Eigen::RowVector4d& p3,
    const Eigen::RowVector4d& p4,
    const Eigen::RowVector4d& p5,
    const double v1,
    const double v2,
    const double v3,
    const double v4,
    const double v5,
    const Eigen::RowVector4d& g1,
    const Eigen::RowVector4d& g2,
    const Eigen::RowVector4d& g3,
    const Eigen::RowVector4d& g4,
    const Eigen::RowVector4d& g5,
    Eigen::RowVector<double, 35>& bezier,
    bool& inside)
{
    // ───── compile–time reciprocals (mul beats div) ───────────────────────────
    constexpr double inv3 = 1.0 / 3.0;
    constexpr double inv6 = 1.0 / 6.0;
    constexpr double half = 0.5;

    // ───── pairwise differences (one evaluation each) ─────────────────────────
    const Eigen::RowVector4d p12 = p1 - p2, p13 = p1 - p3, p14 = p1 - p4, p15 = p1 - p5,
                             p23 = p2 - p3, p24 = p2 - p4, p25 = p2 - p5, p34 = p3 - p4,
                             p35 = p3 - p5, p45 = p4 - p5;

    // ───── per-vertex edge control points (dot once, inv3 mul) ────────────────
    const double v1s0 = v1 - inv3 * g1.dot(p12);
    const double v1s1 = v1 - inv3 * g1.dot(p13);
    const double v1s2 = v1 - inv3 * g1.dot(p14);
    const double v1s3 = v1 - inv3 * g1.dot(p15);

    const double v2s0 = v2 + inv3 * g2.dot(p12);
    const double v2s1 = v2 - inv3 * g2.dot(p23);
    const double v2s2 = v2 - inv3 * g2.dot(p24);
    const double v2s3 = v2 - inv3 * g2.dot(p25);

    const double v3s0 = v3 + inv3 * g3.dot(p13);
    const double v3s1 = v3 + inv3 * g3.dot(p23);
    const double v3s2 = v3 - inv3 * g3.dot(p34);
    const double v3s3 = v3 - inv3 * g3.dot(p35);

    const double v4s0 = v4 + inv3 * g4.dot(p14);
    const double v4s1 = v4 + inv3 * g4.dot(p24);
    const double v4s2 = v4 + inv3 * g4.dot(p34);
    const double v4s3 = v4 - inv3 * g4.dot(p45);

    const double v5s0 = v5 + inv3 * g5.dot(p15);
    const double v5s1 = v5 + inv3 * g5.dot(p25);
    const double v5s2 = v5 + inv3 * g5.dot(p35);
    const double v5s3 = v5 + inv3 * g5.dot(p45);

    // helper macro for tiny triple sums
#define SUM3(a, b, c) ((a) + (b) + (c))

    // ───── face control points (reuse edge sums, no extra temporaries) ────────
    const double e1 = (v1s0 + v1s1 + v2s0 + v2s1 + v3s0 + v3s1) * inv6;
    const double face1 = e1 + (e1 - SUM3(v1, v2, v3) * inv3) * half;

    const double e2 = (v1s0 + v1s2 + v2s0 + v2s2 + v4s0 + v4s1) * inv6;
    const double face2 = e2 + (e2 - SUM3(v1, v2, v4) * inv3) * half;

    const double e3 = (v1s0 + v1s3 + v2s0 + v2s3 + v5s0 + v5s1) * inv6;
    const double face3 = e3 + (e3 - SUM3(v1, v2, v5) * inv3) * half;

    const double e4 = (v1s1 + v1s2 + v3s0 + v3s2 + v4s0 + v4s2) * inv6;
    const double face4 = e4 + (e4 - SUM3(v1, v3, v4) * inv3) * half;

    const double e5 = (v1s1 + v1s3 + v3s0 + v3s3 + v5s0 + v5s2) * inv6;
    const double face5 = e5 + (e5 - SUM3(v1, v3, v5) * inv3) * half;

    const double e6 = (v1s2 + v1s3 + v4s0 + v4s3 + v5s0 + v5s3) * inv6;
    const double face6 = e6 + (e6 - SUM3(v1, v4, v5) * inv3) * half;

    const double e7 = (v2s1 + v2s2 + v3s1 + v3s2 + v4s1 + v4s2) * inv6;
    const double face7 = e7 + (e7 - SUM3(v2, v3, v4) * inv3) * half;

    const double e8 = (v2s1 + v2s3 + v3s1 + v3s3 + v5s1 + v5s2) * inv6;
    const double face8 = e8 + (e8 - SUM3(v2, v3, v5) * inv3) * half;

    const double e9 = (v2s2 + v2s3 + v4s1 + v4s3 + v5s1 + v5s3) * inv6;
    const double face9 = e9 + (e9 - SUM3(v2, v4, v5) * inv3) * half;

    const double e10 = (v3s2 + v3s3 + v4s2 + v4s3 + v5s2 + v5s3) * inv6;
    const double face10 = e10 + (e10 - SUM3(v3, v4, v5) * inv3) * half;

#undef SUM3

    // ───── populate the 35-coeff vector (single << stream) ────────────────────
    bezier << v1, v2, v3, v4, v5, v1s0, v1s1, v1s2, v1s3, v2s0, v2s1, v2s2, v2s3, v3s0, v3s1, v3s2,
        v3s3, v4s0, v4s1, v4s2, v4s3, v5s0, v5s1, v5s2, v5s3, face1, face2, face3, face4, face5,
        face6, face7, face8, face9, face10;

#ifndef no_inside
    inside =
        inside || bezier({0, 1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 25, 26, 28, 31})
                          .maxCoeff() <= 0;
    inside = inside ||
             bezier({1, 2, 3, 4, 10, 11, 12, 14, 15, 16, 18, 19, 20, 22, 23, 24, 31, 32, 33, 34})
                     .maxCoeff() <= 0;
#endif

    // ───── final sign-test (unchanged logic) ──────────────────────────────────
    const double mx = bezier.maxCoeff();
    const double mn = bezier.minCoeff();
    return get_sign(mx) != get_sign(mn);
}

/// Compute the adjugate (classical adjoint) of a 4×4 matrix using explicit minor and cofactor expansion.
/// The adjugate matrix is the transpose of the cofactor matrix and is used in computing the matrix inverse
/// or determinant-related geometric transformations in 4D computations. This implementation iteratively
/// evaluates all 3×3 minors of the input matrix and constructs the transposed cofactor matrix explicitly.
///
/// @param[in] mat          The input 4×4 matrix. Typically represents a Jacobian or transformation matrix
///                         in 4D space (e.g., spatial-temporal derivatives or local frame transformations).
/// @param[out] adjugate    The resulting 4×4 adjugate (classical adjoint) matrix. Each element (j, i)
///                         corresponds to the signed determinant of the 3×3 minor obtained by removing
///                         the i-th row and j-th column from the input matrix.
void adjugate(const Eigen::Matrix<double, 4, 4>& mat, Eigen::Matrix4d& adjugate)
{
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            Eigen::Matrix3d minor;

            int rowOffset = 0;
            for (int row = 0; row < 4; ++row) {
                if (row == i) {
                    rowOffset = 1;
                    continue;
                }
                int colOffset = 0;
                for (int col = 0; col < 4; ++col) {
                    if (col == j) {
                        colOffset = 1;
                        continue;
                    }
                    minor(row - rowOffset, col - colOffset) = mat(row, col);
                }
            }
            // Compute the cofactor
            double cofactor = std::pow(-1, i + j) * minor.determinant();
            adjugate(j, i) = cofactor; // Note the transpose here
        }
    }
}

/// Perform degree elevation from a quadratic Bezier simplex to a cubic Bezier simplex in 4D space.
/// This function computes the elevated Bezier control ordinates by linearly combining the 16 original
/// cubic coefficients using a precomputed elevation matrix (`elevMatrix`) and weighting matrix (`ls`).
/// The resulting 35 cubic control ordinates preserve the original geometric form while improving
/// smoothness and compatibility with higher-order Bezier simplex evaluations.
///
/// @param[in] ords         The 16 Bezier control ordinates (coefficients) corresponding to the cubic Bezier tetrahedron.
///                         Each ordinate represents the scalar value of the Bezier function at a vertex in barycentric order.
/// @param[out] bezier      The 35 elevated Bezier control ordinates for the cubic Bezier simplex.
///                         These ordinates are computed by applying the linear transformation defined by `elevMatrix`
///                         and `ls`, followed by normalization.
void bezierElev(const Eigen::RowVector<double, 16>& ords, Eigen::RowVector<double, 35>& bezier)
{
    Eigen::Vector<double, 5> elevRow;
    for (size_t i = 0; i < 35; i++) {
        elevRow << ords[elevMatrix[i][0]], ords[elevMatrix[i][1]], ords[elevMatrix[i][2]],
            ords[elevMatrix[i][3]], ords[elevMatrix[i][4]];
        bezier[i] = ls.row(i).dot(elevRow);
    }
    bezier = bezier / 3.0;
}

/// Compute the elevated Bezier control ordinates for the directional derivative of a cubic (order-3)
/// 4D Bezier simplex along its temporal edge. The derivative yields a quadratic (order-2) Bezier simplex,
/// which is subsequently degree-elevated to cubic to produce 35 ordinates for uniform processing.
///
/// @param[in] ords         35 cubic Bezier ordinates of the scalar field on the 4D 5-cell (barycentric order).
/// @param[in] verts        The 5 vertices of the simplex (RowVector4d), used to normalize by ‖verts[0]−verts[4]‖.
/// @param[out] bezierGrad  35 cubic Bezier ordinates of the directional derivative after degree elevation.
///
/// @details                Differences indexed by `derMatrix` form 16 cubic derivative ordinates (scaled by 3/‖Δt‖),
///                         then `bezierElev` elevates quadratic→cubic. Sign consistency of `bezierGrad` is returned.
///
/// @return                 true if max and min coefficients of `bezierGrad` share the same sign (no zero crossing);
///                         false otherwise.
bool bezierDerOrds(
    const Eigen::RowVector<double, 35>& ords,
    const std::array<Eigen::RowVector4d, 5>& verts,
    Eigen::RowVector<double, 35>& bezierGrad)
{
    Eigen::RowVector<double, 16> vals;
    double norm = (verts[0] - verts[4]).norm();
    // Loop through rows of derMatrix
    for (int i = 0; i < 15; i++) {
        vals[i] = (ords[derMatrix[i][0]] - ords[derMatrix[i][1]]) * 3.0 /
                  norm; // Normalize by the norm of verts
    }
    vals[15] = 0;
    bezierElev(vals, bezierGrad);
    return get_sign(bezierGrad.maxCoeff()) == get_sign(bezierGrad.minCoeff());
}

/// Extremely fast outside-hull clipping test for a 2D set of Bezier ordinates projected to the plane.
/// Given 35 2D points (e.g., a cubic Bezier simplex control points in 4D),
/// this routine checks whether all points lie within some origin-centered half-plane
/// (i.e., there exists a direction u such that u·p >= 0 for all points p).
/// Equivalently, it verifies that the convex hull of the points does **not** wrap around
/// the origin (no full 360° coverage), so a line through the origin separates the origin
/// from the point set.
///
/// The algorithm maintains a feasible angular wedge of outward normals, initialized by the
/// first point’s perpendiculars, and iteratively clips this wedge with each subsequent point.
/// If the wedge becomes empty, the origin is enclosed by the hull and the test fails.
///
/// @param[in] pts          A 2×35 matrix whose columns are the 2D points to test.
///                         Columns are 35 bezier control points.
/// @return                 true  if a separating origin-centered half-plane exists (origin is outside the hull);
///                         false if the feasible wedge collapses (origin is enclosed/covered by the hull).
///
/// @note                   Robustness control: uses a small epsilon (1e-7) to handle near-collinear cases and ties by
///                         falling back to an orientation test via a perpendicular vector.
///                         The helper `perp(v) = (-v_y, v_x)` rotates vectors by +90°.
bool outHullClip2D(Eigen::Matrix<double, 2, 35> pts)
{
    constexpr double eps = 0.0000001;
    bool r1, r2;
    double t;
    auto perp = [](const Eigen::Vector2d& data) { return Eigen::Vector2d{-data[1], data[0]}; };
    std::array<Eigen::Vector2d, 2> range = {-perp(pts.col(0)), perp(pts.col(0))};
    for (size_t i = 1; i < 35; i++) {
        t = range[0].dot(pts.col(i));
        if (t > eps) {
            r1 = true;
        } else if (t < -eps) {
            r1 = false;
        } else if (perp(range[0]).dot(pts.col(i)) > 0) {
            r1 = true;
        } else {
            return false;
        }
        t = range[1].dot(pts.col(i));
        if (t > eps) {
            r2 = true;
        } else if (t < -eps) {
            r2 = false;
        } else if (perp(range[1]).dot(pts.col(i)) < 0) {
            r2 = true;
        } else {
            return false;
        }

        if (!r1 && !r2) {
            return false;
        } else if (!r1) {
            range[0] = -perp(pts.col(i));
        } else if (!r2) {
            range[1] = perp(pts.col(i));
        }
    }
    return true;
}

bool refineFtBezier(
    const std::array<vertex4d*, 5>& verts,
    const double threshold,
    bool& inside,
    bool& choice,
    bool& zeroX)
{
    const auto& p1 = verts[0]->coord;
    const auto& p2 = verts[1]->coord;
    const auto& p3 = verts[2]->coord;
    const auto& p4 = verts[3]->coord;
    const auto& p5 = verts[4]->coord;

    const auto& v1 = verts[0]->valGradList.first;
    const auto& v2 = verts[1]->valGradList.first;
    const auto& v3 = verts[2]->valGradList.first;
    const auto& v4 = verts[3]->valGradList.first;
    const auto& v5 = verts[4]->valGradList.first;

    const auto& g1 = verts[0]->valGradList.second;
    const auto& g2 = verts[1]->valGradList.second;
    const auto& g3 = verts[2]->valGradList.second;
    const auto& g4 = verts[3]->valGradList.second;
    const auto& g5 = verts[4]->valGradList.second;
    Eigen::RowVector<double, 35> bezierVals;
    if (!bezier4D(p1, p2, p3, p4, p5, v1, v2, v3, v4, v5, g1, g2, g3, g4, g5, bezierVals, inside)) {
        return false;
    };
    return true;
}

/// See header
bool refineFt(
    const std::array<vertex4d*, 5>& verts,
    const double threshold,
    bool& inside,
    bool& choice,
    bool& zeroX,
    std::array<double, timer_amount>& profileTimer,
    std::array<size_t, timer_amount>& profileCount)
{
#if time_profile
    Timer first_bezier_timer(first_bezier, [&](auto timer, auto ms) {
        combine_timer(profileTimer, profileCount, timer, ms);
    });
#endif
    const auto& p1 = verts[0]->coord;
    const auto& p2 = verts[1]->coord;
    const auto& p3 = verts[2]->coord;
    const auto& p4 = verts[3]->coord;
    const auto& p5 = verts[4]->coord;

    const auto& v1 = verts[0]->valGradList.first;
    const auto& v2 = verts[1]->valGradList.first;
    const auto& v3 = verts[2]->valGradList.first;
    const auto& v4 = verts[3]->valGradList.first;
    const auto& v5 = verts[4]->valGradList.first;

    const auto& g1 = verts[0]->valGradList.second;
    const auto& g2 = verts[1]->valGradList.second;
    const auto& g3 = verts[2]->valGradList.second;
    const auto& g4 = verts[3]->valGradList.second;
    const auto& g5 = verts[4]->valGradList.second;
    Eigen::RowVector<double, 35> bezierVals;
    if (!bezier4D(p1, p2, p3, p4, p5, v1, v2, v3, v4, v5, g1, g2, g3, g4, g5, bezierVals, inside)) {
#if time_profile
        first_bezier_timer.Stop();
#endif
        return false;
    };
    Eigen::RowVector<double, 35> bezierGrad;
    if (bezierDerOrds(bezierVals, {p1, p2, p3, p4, p5}, bezierGrad)) {
#if time_profile
        first_bezier_timer.Stop();
#endif
        return false;
    }
    Eigen::Matrix<double, 2, 35> nPoints_eigen;
    nPoints_eigen << bezierVals, bezierGrad;
    zeroX = !outHullClip2D(nPoints_eigen);
#if time_profile
    first_bezier_timer.Stop();
    Timer first_func_timer(first_func, [&](auto timer, auto ms) {
        combine_timer(profileTimer, profileCount, timer, ms);
    });
#endif
    if (zeroX) {
        auto vec1 = p2 - p1, vec2 = p3 - p1, vec3 = p4 - p1, vec4 = p5 - p1;
        Eigen::Matrix4d vec;
        vec << vec1, vec2, vec3, vec4;
        Eigen::Matrix4d adj;
        adjugate(vec, adj);
        auto v1 = bezierGrad[0];
        auto v2 = bezierGrad[1];
        auto v3 = bezierGrad[2];
        auto v4 = bezierGrad[3];
        auto v5 = bezierGrad[4];
        Eigen::RowVector4d gradList = Eigen::RowVector4d(v2 - v1, v3 - v1, v4 - v1, v5 - v1) * adj;
        Eigen::RowVector<double, 30> error =
            ((bezierGrad.tail(30) - (bezierGrad.head(5) * ls.bottomRows(30).transpose()) / 3) *
             vec.determinant())
                .array()
                .abs();
        if (error.maxCoeff() > threshold * gradList.norm()) {
            Eigen::RowVector<double, 16> topFError = error(topFIndices);
            Eigen::RowVector<double, 16> botFError = error(botFIndices);
            choice = std::max(error[3], error[16]) >
                     std::min(topFError.maxCoeff(), botFError.maxCoeff());
#if time_profile
            first_func_timer.Stop();
#endif
            return true;
        }
    }
#if time_profile
    first_func_timer.Stop();
#endif
    return false;
}

/// Construct the values of one function at the bezier control points within a tet.
///
/// @param[in] p1-p4          The vertex cooridantes of four tet vertices.
/// @param[in] v1-v4         The function values at four tet vertices.
/// @param[in] g1-g4            The total derivative of the functions in x, y, z, t direction at four tet vertices.
/// @param[out] bezier          The eigen vector of 20 bezier values.
///
/// @return         If 20 bezier values contain zero-crossing.
bool bezier3D(
    const Eigen::RowVector4d& p1,
    const Eigen::RowVector4d& p2,
    const Eigen::RowVector4d& p3,
    const Eigen::RowVector4d& p4,
    const double v1,
    const double v2,
    const double v3,
    const double v4,
    const Eigen::RowVector4d& g1,
    const Eigen::RowVector4d& g2,
    const Eigen::RowVector4d& g3,
    const Eigen::RowVector4d& g4,
    Eigen::RowVector<double, 20>& bezier)
{
    // Compute edge values
    std::array<double, 4> v1s = {
        v1 + g1.dot(p2 - p1) / 3,
        v1 + g1.dot(p3 - p1) / 3,
        v1 + g1.dot(p4 - p1) / 3};

    std::array<double, 4> v2s = {
        v2 + g2.dot(p3 - p2) / 3,
        v2 + g2.dot(p4 - p2) / 3,
        v2 + g2.dot(p1 - p2) / 3};

    std::array<double, 4> v3s = {
        v3 + g3.dot(p4 - p3) / 3,
        v3 + g3.dot(p1 - p3) / 3,
        v3 + g3.dot(p2 - p3) / 3};

    std::array<double, 4> v4s = {
        v4 + g4.dot(p1 - p4) / 3,
        v4 + g4.dot(p2 - p4) / 3,
        v4 + g4.dot(p3 - p4) / 3};
    // Compute face values
    double face1 =
        (9 * (v2s[0] + v2s[1] + v3s[0] + v3s[2] + v4s[1] + v4s[2]) / 6 - v2 - v3 - v4) / 6;
    double face2 =
        (9 * (v1s[1] + v1s[2] + v3s[0] + v3s[1] + v4s[0] + v4s[2]) / 6 - v1 - v3 - v4) / 6;
    double face3 =
        (9 * (v1s[0] + v1s[2] + v2s[1] + v2s[2] + v4s[0] + v4s[1]) / 6 - v1 - v2 - v4) / 6;
    double face4 =
        (9 * (v1s[0] + v1s[1] + v2s[0] + v2s[2] + v3s[1] + v3s[2]) / 6 - v1 - v2 - v3) / 6;

    // Combine results into a single row vector
    bezier << v1, v2, v3, v4, v1s[0], v1s[1], v1s[2], v2s[0], v2s[1], v2s[2], v3s[0], v3s[1],
        v3s[2], v4s[0], v4s[1], v4s[2], face1, face2, face3, face4;

    if (get_sign(bezier.maxCoeff()) == get_sign(bezier.minCoeff())) {
        return false;
    }
    return true;
}

/// Construct the value differences between linear interpolations and bezier approximations at 16 bezier control points (excluding control points at tet vertices)
/// @param[in] valList          The eigen vector of 20 bezier values.
///
/// @return         The value differences at 16 control points.
Eigen::Vector<double, 16> bezierDiff(const Eigen::Vector<double, 20> valList)
{
    /// Constant coefficient to obtain linear interpolated values at each bezier control points
    const Eigen::Matrix<double, 16, 4> linear_coeff{
        {2, 1, 0, 0},
        {2, 0, 1, 0},
        {2, 0, 0, 1},
        {0, 2, 1, 0},
        {0, 2, 0, 1},
        {1, 2, 0, 0},
        {0, 0, 2, 1},
        {1, 0, 2, 0},
        {0, 1, 2, 0},
        {1, 0, 0, 2},
        {0, 1, 0, 2},
        {0, 0, 1, 2},
        {0, 1, 1, 1},
        {1, 0, 1, 1},
        {1, 1, 0, 1},
        {1, 1, 1, 0}};
    Eigen::Vector<double, 16> linear_val = (linear_coeff * valList.head(4)) / 3;
    return valList.tail(16) - linear_val;
}

/// See header
bool refine3D(const std::array<vertex4d, 4>& verts, const double threshold)
{
    const auto& p1 = verts[0].coord;
    const auto& p2 = verts[1].coord;
    const auto& p3 = verts[2].coord;
    const auto& p4 = verts[3].coord;

    const auto& v1 = verts[0].valGradList.first;
    const auto& v2 = verts[1].valGradList.first;
    const auto& v3 = verts[2].valGradList.first;
    const auto& v4 = verts[3].valGradList.first;

    const auto& g1 = verts[0].valGradList.second;
    const auto& g2 = verts[1].valGradList.second;
    const auto& g3 = verts[2].valGradList.second;
    const auto& g4 = verts[3].valGradList.second;

    Eigen::RowVector<double, 20> bezierVals;
    if (!bezier3D(p1, p2, p3, p4, v1, v2, v3, v4, g1, g2, g3, g4, bezierVals)) {
        return false;
    };
    Eigen::Vector3d vec1 = p2.head(3) - p1.head(3), vec2 = p3.head(3) - p1.head(3),
                    vec3 = p4.head(3) - p1.head(3);
    Eigen::Matrix3d vec;
    vec << vec1, vec2, vec3;
    double D = vec.determinant();
    Eigen::Matrix3d crossMatrix;
    crossMatrix << vec2.cross(vec3), vec3.cross(vec1), vec1.cross(vec2);
    Eigen::Vector3d unNormF =
        Eigen::RowVector3d(v2 - v1, v3 - v1, v4 - v1) * crossMatrix.transpose();
    double lhs = (bezierDiff(bezierVals) * D).array().abs().maxCoeff();
    double rhs = threshold * unNormF.norm();
    if (lhs > rhs) {
        return true;
    }
    return false;
}
