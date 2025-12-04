#include <span>
#include <queue>
#include <optional>
#include <CLI/CLI.hpp>
#include <chrono>
#include <igl/read_triangle_mesh.h>
#include <igl/signed_distance.h>
#include <nlohmann/json.hpp>
#include <algorithm>

#include <lagrange/io/save_mesh.h>
#include <lagrange/mesh_cleanup/remove_degenerate_facets.h>
#include <lagrange/mesh_cleanup/remove_topologically_degenerate_facets.h>
#include <lagrange/views.h>
#include <lagrange/topology.h>
#include <lagrange/utils/SmallVector.h>

#include <generalized_sweep/generalized_sweep.h>


//#include "init_grid.h"
//#include "io.h"
//#include "col_gridgen.h"
#include "trajectory.h"
//#include "post_processing.h"
//#include "timer.h"

#define SAVE_CONTOUR 0
#define batch_stats 0
#define batch_time 0

template <typename Scalar, typename Index>
void save_features(std::string_view filename, lagrange::SurfaceMesh<Scalar, Index>& arrangement)
{
    std::ofstream fout(filename.data());
    if (!fout) {
        throw std::runtime_error("Failed to open file: " + std::string(filename));
    }

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

int main(int argc, const char *argv[])
{
    struct
    {
        std::string grid_file;
        std::string output_path;
        std::string function_file = "";
        double threshold = 0.0005;
        double traj_threshold = 0.005;
        int max_splits = std::numeric_limits<int>::max();
        int rot = 0;
        bool without_insideness_check = false;
        bool without_snapping = false;
        bool without_opt_triangulation = false;
        bool cyclic = false;
    } args;
    CLI::App app{"Generalized Swept Volume"};
    app.add_option("grid", args.grid_file, "Initial grid file")->required();
    app.add_option("output", args.output_path, "Output path")->required();
    app.add_option("-f,--function", args.function_file, "Implicit function file");
    app.add_option("--ee,--epsilon-env", args.threshold, "Envelope threshold");
    app.add_option("--es, --epsilon-sil", args.traj_threshold, "Silhouette threshold");
    app.add_flag("--without-inside-check", args.without_insideness_check, "Turn on the refinement for the inside regions of the envelope");
    app.add_option("-m,--max-splits", args.max_splits, "Maximum number of splits");
    app.add_option("-r,--rotation-number", args.rot, "Number of rotations");
    app.add_flag("--without-snapping", args.without_snapping, "Disable vertex snapping in iso-surfacing step");
    app.add_flag("--without-optimal-triangulation", args.without_opt_triangulation,
                 "Disable optimal triangulation in iso-surfacing triangulation step");
    app.add_flag("--cyclic", args.cyclic, "Whether the trajectory is cyclic or not");
    CLI11_PARSE(app, argc, argv);

    using Scalar = sweep::Scalar;


    std::string output_path = args.output_path;
    int max_splits = args.max_splits;
    bool insideness_check = !args.without_insideness_check;
    std::string function_file = args.function_file;
    double threshold = args.threshold;
    double traj_threshold = args.traj_threshold;
    int rotation = args.rot;
    const int dim = 4;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::AABB<Eigen::MatrixXd,3> tree;
    igl::FastWindingNumberBVH fwn_bvh;
    std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> implicit_sweep;
    if (args.function_file != ""){
        if (args.function_file.find(".obj") != std::string::npos){
            igl::read_triangle_mesh(args.function_file,V,F);
            tree.init(V,F);
            int order = 2;
            igl::fast_winding_number(V, F, order, fwn_bvh);\
            igl::WindingNumberAABB<double, int> hier;
            hier.set_mesh(V,F);
            hier.grow();
            /// main function:
            /// the lambda function for function evaluations
            ///  @param[in] data            The 4D coordinate
            ///  @return    A std::pari<Scalar, Eigen::RowVector4d> of the value and the gradients at this 4D point
            ///
            ///libigl input using mesh files (unstable gradients, need high resolution mesh input):
            implicit_sweep = [&](Eigen::RowVector4d data)->std::pair<Scalar, Eigen::RowVector4d>
            {
                Scalar value;
                Eigen::RowVector4d gradient;
                const double iso = 0.001;
                Eigen::RowVector3d P = data.head(3);
                double t = data[3];
                Eigen::RowVector3d running_closest_point = V.row(0);
                double running_sign = 1.0;
                int i;
                double s,sqrd,sqrd2,s2;
                Eigen::Matrix3d VRt,Rt;
                Eigen::RowVector3d xt,vt,pos,c,c2,xyz_grad,point_velocity;
                trajLine3D2(t, xt, vt);
                trajLineRot3D(t, Rt, VRt, rotation);
                pos = ((Rt.inverse())*((P - xt).transpose())).transpose();
                // fast winding number
                Eigen::VectorXd w;
                igl::fast_winding_number(fwn_bvh,2.0,pos,w);
                s = 1.-2.*w(0);
                sqrd = tree.squared_distance(V,F,pos,i,c);
                value = s*sqrt(sqrd);
                Eigen::RowVector3d cp = c - pos;
                cp.normalize();
                xyz_grad  = (-s) * cp * Rt.inverse();
                gradient.template head<3>() << xyz_grad;
                point_velocity = (-Rt.inverse()*VRt*Rt.inverse()*(P.transpose() - xt.transpose()) - Rt.inverse()*vt.transpose()).transpose();
                gradient(3) =  (-s) * cp.dot(point_velocity);
                return {value, gradient};
            };
        } else if (args.function_file == "elbow") {
            implicit_sweep = elbow;
        } else if (args.function_file == "bezier") {
            implicit_sweep = bezier;
        } else if (args.function_file == "blend_spheres") {
            implicit_sweep = blend_spheres;
        } else if (args.function_file == "blend_sphere_torus") {
            implicit_sweep = blend_sphere_torus;
        } else if (args.function_file == "sphere_spiral") {
            implicit_sweep = sphere_spiral;
        } else if (args.function_file == "knot") {
            implicit_sweep = knot;
        } else if (args.function_file == "brush_stroke") {
            implicit_sweep = brush_stroke;
        } else if (args.function_file == "brush_stroke_blending") {
            implicit_sweep = brush_stroke_blending;
        } else if (args.function_file == "concentric_rings") {
            implicit_sweep = concentric_rings;
        } else if (args.function_file == "spinning_rod") {
            implicit_sweep = spinning_rod;
        } else if (args.function_file == "letter_L") {
            implicit_sweep = letter_L;
        } else if (args.function_file == "letter_L_blend") {
            implicit_sweep = letter_L_blend;
        } else if (args.function_file == "torus_rotation") {
            implicit_sweep = torus_rotation;
        } else if (args.function_file == "loopDloop_with_offset") {
            implicit_sweep = loopDloop_with_offset;
        } else if (args.function_file == "loopDloop_with_offset_v2") {
            implicit_sweep = loopDloop_with_offset_v2;
        } else if (args.function_file == "doghead") {
            implicit_sweep = doghead;
        } else if (args.function_file == "star_S") {
            implicit_sweep = star_S;
        } else if (args.function_file == "star_D") {
            implicit_sweep = star_D;
        } else if (args.function_file == "star_F") {
            implicit_sweep = star_F;
        } else if (args.function_file == "star_I") {
            implicit_sweep = star_I;
        } else if (args.function_file == "fertility") {
            implicit_sweep = fertility;
        } else if (args.function_file == "fertility_v2") {
            implicit_sweep = fertility_v2;
        } else if (args.function_file == "fertility_v3") {
            implicit_sweep = fertility_v3;
        } else if (args.function_file == "fertility_v4") {
            implicit_sweep = fertility_v4;
        } else if (args.function_file == "fertility_v5") {
            implicit_sweep = fertility_v5;
        } else if (args.function_file == "fertility_v6") {
            implicit_sweep = fertility_v6;
        } else if (args.function_file == "bunny_blend") {
            implicit_sweep = bunny_blend;
        } else if (args.function_file == "loopDloop_with_offset_v3") {
            implicit_sweep = loopDloop_with_offset_v3;
        } else if (args.function_file == "VIPSS_blend") {
            implicit_sweep = VIPSS_blend;
        } else if (args.function_file == "VIPSS_blend3") {
            implicit_sweep = VIPSS_blend3;
        } else if (args.function_file == "VIPSS_blend_S") {
            implicit_sweep = VIPSS_blend_S;
        } else if (args.function_file == "wheel_I") {
            implicit_sweep = wheel_I;
        } else if (args.function_file == "wheel_I_shrink") {
            implicit_sweep = wheel_I_shrink;
        } else if (args.function_file == "mesh_I") {
            implicit_sweep = mesh_I;
        } else if (args.function_file == "tangle_cube_roll") {
            implicit_sweep = tangle_cube_roll;
        } else if (args.function_file == "ball_genus_roll") {
            implicit_sweep = ball_genus_roll;
        } else if (args.function_file == "tangle_chair_S") {
            implicit_sweep = tangle_chair_S;
        } else if (args.function_file == "tangle_chair") {
            implicit_sweep = tangle_chair;
        } else if (args.function_file == "rotating_rod") {
            implicit_sweep = rotating_rod;
        } else if (args.function_file == "nested_balls") {
            implicit_sweep = nested_balls;
        } else if (args.function_file == "close_loop") {
            implicit_sweep = close_loop;
        } else if (args.function_file == "close_loop_2") {
            implicit_sweep = close_loop_2;
        } else if (args.function_file == "tet_roll") {
            implicit_sweep = tet_roll;
        } else {
            throw std::runtime_error("ERROR: file format not supported");
        }
    }else{
        /// use hard coded models as default/testing purpose.
        implicit_sweep = [&](Eigen:: RowVector4d data){
            return flippingDonutFullTurn(data);
        };
    }


    // TODO: refactor into a function
    sweep::GridSpec grid_spec;
    {
        using json = nlohmann::json;
        std::ifstream fin(args.grid_file.c_str());
        if (!fin) {
            throw std::runtime_error("tet mesh file not exist!");
        }
        json data;
        fin >> data;
        fin.close();
        if (!data.contains("resolution") || !data.contains("bbox_min") || !data.contains("bbox_max")) {
            throw std::runtime_error("grid specification missing in json file!");
        }
        if (data["resolution"].size() == 3) {
            grid_spec.resolution = {
                data["resolution"][0].get<size_t>(),
                data["resolution"][1].get<size_t>(),
                data["resolution"][2].get<size_t>()
            };
        } else {
            assert(data["resolution"].size() == 1);
            size_t res = data["resolution"][0].get<size_t>();
            grid_spec.resolution = {res, res, res};
        }
        grid_spec.bbox_min = {
            data["bbox_min"][0].get<float>(),
            data["bbox_min"][1].get<float>(),
            data["bbox_min"][2].get<float>()
        };
        grid_spec.bbox_max = {
            data["bbox_max"][0].get<float>(),
            data["bbox_max"][1].get<float>(),
            data["bbox_max"][2].get<float>()
        };
    }

    sweep::SweepOptions options;
    options.epsilon_env = threshold;
    options.epsilon_sil = traj_threshold;
    options.max_split = max_splits;
    options.with_insideness_check = insideness_check;
    options.with_snapping = !args.without_snapping;
    options.cyclic = args.cyclic;
    auto result = sweep::generalized_sweep(implicit_sweep, std::move(grid_spec),
            std::move(options));
    auto& envelope = result.envelope;
    auto& sweep_surface = result.sweep_surface;
    auto& sweep_arrangement = result.arrangement;


    // Saving result
    auto saving_start = std::chrono::time_point_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now()
    ).time_since_epoch().count();

    if (!std::filesystem::exists(output_path)) {
        // Attempt to create the directory
        if (std::filesystem::create_directory(output_path)) {
            std::cout << "Directory created successfully." << std::endl;
        } else {
            std::cerr << "Failed to create directory." << std::endl;
        }
    } else {
        std::cout << "Directory already exists. Removing all contents..." << std::endl;
        for (const auto& entry : std::filesystem::directory_iterator(output_path)) {
            std::error_code ec;
            std::filesystem::remove_all(entry.path(), ec);
            if (ec) {
                std::cerr << "Error removing " << entry.path() << ": " << ec.message() << std::endl;
            }
        }
    }

    lagrange::io::save_mesh(output_path + "/envelope.msh", envelope);
    lagrange::io::save_mesh(output_path + "/sweep_surface.msh", sweep_surface);
    lagrange::io::save_mesh(output_path + "/arrangement.msh", sweep_arrangement);
    save_features(output_path + "/features.obj", sweep_arrangement);
#if SAVE_CONTOUR
    //mtet::save_mesh(output_path + "/tet_grid.msh", grid);
    //save_grid_for_mathematica(output_path + "/contour_iso.json", grid, vertexMap);
#endif

    auto saving_end = std::chrono::time_point_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now()
    ).time_since_epoch().count();
    std::cout << "Saving time: "
        << (saving_end - saving_start) * 1e-6 << " seconds" << std::endl;

    return 0;
}
