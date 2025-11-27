#include <span>
#include <queue>
#include <optional>
#include <CLI/CLI.hpp>
#include <mtet/io.h>
#include <chrono>
#include <igl/read_triangle_mesh.h>
#include <igl/signed_distance.h>
#include <mtetcol/contour.h>
#include <mtetcol/simplicial_column.h>
#include <mtetcol/io.h>
#include <igl/write_triangle_mesh.h>
#include <igl/remove_unreferenced.h>
#include <nlohmann/json.hpp>
#include <algorithm>

#include <lagrange/io/save_mesh.h>
#include <lagrange/mesh_cleanup/remove_degenerate_facets.h>
#include <lagrange/mesh_cleanup/remove_topologically_degenerate_facets.h>
#include <lagrange/views.h>
#include <lagrange/topology.h>
#include <lagrange/utils/SmallVector.h>
#include <lagrange/triangulate_polygonal_facets.h>


#include "init_grid.h"
#include "io.h"
#include "col_gridgen.h"
#include "trajectory.h"
#include "post_processing.h"
#include "timer.h"

#define SAVE_CONTOUR 1
#define batch_stats 0
#define batch_time 0

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
    // Read initial grid
    mtet::MTetMesh grid;
    if (args.grid_file.find(".json") != std::string::npos){
        grid = init_grid::load_tet_mesh(args.grid_file);
        mtet::save_mesh("init.msh", grid);
        grid = mtet::load_mesh("init.msh");
    } else {
        grid = mtet::load_mesh(args.grid_file);
    }
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
    ///
    ///
    ///Grid generation:
    vertExtrude vertexMap;
    insidenessMap insideMap;
    std::cout << "Start to generate the background grid..." << std::endl;
    std::array<double, timer_amount> profileTimer{};
    std::array<size_t, timer_amount> profileCount{};
    auto starterTime = std::chrono::high_resolution_clock::now();
    spdlog::set_level(spdlog::level::off);
    if (!gridRefine(grid, vertexMap, insideMap, implicit_sweep, threshold, traj_threshold, max_splits, insideness_check, profileTimer, profileCount)){
        throw std::runtime_error("ERROR: grid generation failed");
        return 0;
    };
    spdlog::set_level(spdlog::level::info);
    
    Scalar iso_value = 0.0;
    bool cyclic = args.cyclic;
    std::vector<mtetcol::Scalar> verts;
    std::vector<mtetcol::Index> simps;
    std::vector<std::vector<double>> time;
    std::vector<std::vector<double>> values;
    convert_4d_grid_mtetcol(grid, vertexMap,
                            verts,
                            simps,
                            time,
                            values,
                            cyclic);
    auto stopperTime = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::time_point_cast<std::chrono::microseconds>(starterTime).time_since_epoch().count();
    auto grid_end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    std::cout << "Grid-building time: " << (grid_end - start) * 1e-6 << " seconds" << std::endl;
    std::function<std::span<double>(size_t)> time_func = [&](size_t index)->std::span<double>{
        return time[index];
    };
    std::function<std::span<double>(size_t)> values_func = [&](size_t index)->std::span<double>{
        return values[index];
    };
    mtetcol::SimplicialColumn<dim> columns;
    columns.set_vertices(verts);
    columns.set_simplices(simps);
    columns.set_time_samples(time_func, values_func);
    
    auto contour = columns.extract_contour(iso_value, cyclic);
    if (!contour.is_manifold()) {
        throw std::runtime_error("ERROR: extracted contour is not manifold");
    }
    stopperTime = std::chrono::high_resolution_clock::now();
    auto surface_1_end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    std::cout << "Surfacing time: " << (surface_1_end - grid_end) * 1e-6 << " seconds (First marching)" << std::endl;
    size_t num_contour_vertices = contour.get_num_vertices();
    std::vector<double> function_values(num_contour_vertices);
    std::vector<double> gradient_values(num_contour_vertices * dim);
    for (uint32_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        auto pos_eval = implicit_sweep(Eigen::RowVector4d{pos[0], pos[1], pos[2], pos[3]});
        function_values[i] = pos_eval.first;
        gradient_values[dim * i] = pos_eval.second[0];
        gradient_values[dim * i + 1] = pos_eval.second[1];
        gradient_values[dim * i + 2] = pos_eval.second[2];
        gradient_values[dim * i + 3] = pos_eval.second[3];
    }
    // Extract isocontour
    auto isocontour = contour.isocontour(function_values, gradient_values, !args.without_snapping);
    if (!isocontour.is_manifold()) {
        std::cout << "isocontour problem" << std::endl;
        throw std::runtime_error("ERROR: extracted isocontour is not manifold");
    }
    //isocontour.triangulate_cycles(!args.without_opt_triangulation);
    //if (!isocontour.is_manifold()) {
    //    std::cout << "triangulated isocontour problem" << std::endl;
    //    throw std::runtime_error("ERROR: extracted isocontour is not manifold");
    //}

    lagrange::SurfaceMesh<double, uint32_t> envelope;
    {
        // Port the isocontour into lagrange mesh
        size_t num_vertices = isocontour.get_num_vertices();
        size_t num_cycles= isocontour.get_num_cycles();

        // Add vertices and time
        envelope.add_vertices(num_vertices);
        envelope.template create_attribute<double>(
            "time",
            lagrange::AttributeElement::Vertex,
            lagrange::AttributeUsage::Scalar,
            1
        );
        auto time_values = attribute_vector_ref<double>(envelope, "time");

        for (size_t i=0; i<num_vertices; i++) {
            auto xyzt = isocontour.get_vertex(i);
            auto pos = envelope.ref_position(i);
            pos[0] = xyzt[0];
            pos[1] = xyzt[1];
            pos[2] = xyzt[2];
            time_values[i] = xyzt[3];
        }

        ankerl::unordered_dense::map<std::pair<mtetcol::Index, mtetcol::Index>, std::vector<size_t>> edge_valence_map;

        // Add polygons
        lagrange::SmallVector<uint32_t, 16> polygon;
        for (size_t i=0; i<num_cycles; i++) {
            auto cycle = isocontour.get_cycle(i);
            size_t cycle_size = cycle.size();
            polygon.clear();
            polygon.resize(cycle_size);

            size_t ind = 0;
            for (auto si : cycle) {
                mtetcol::Index seg_id = index(si);
                bool seg_ori = mtetcol::orientation(si);
                auto seg = isocontour.get_segment(seg_id);
                polygon[ind] = (seg_ori ? seg[0] : seg[1]);
                std::pair<mtetcol::Index, mtetcol::Index> edge_key = {
                    std::min(seg[0], seg[1]),
                    std::max(seg[0], seg[1])
                };

                if (edge_valence_map.find(edge_key) == edge_valence_map.end()) {
                    edge_valence_map[edge_key] = {};
                }
                edge_valence_map[edge_key].push_back(i);

                ind ++;
            }
            envelope.add_polygon({polygon.data(), polygon.size()});
        }

        for (const auto& [edge, cycle_ids] : edge_valence_map) {
            if (cycle_ids.size() > 2) {
                std::cout << "Non-manifold edge detected in the envelope mesh." << std::endl;
                std::cout << "Edge " << edge.first << ", " << edge.second << " is shared by "
                          << cycle_ids.size() << " cycles:" << std::endl;
                for (auto cid : cycle_ids) {
                    std::cout << "  Cycle " << cid << ": ";
                    auto cycle = isocontour.get_cycle(cid);
                    for (auto si : cycle) {
                        mtetcol::Index seg_id = index(si);
                        bool seg_ori = mtetcol::orientation(si);
                        auto seg = isocontour.get_segment(seg_id);
                        std::cout << (seg_ori ? "+" : "-") << seg_id + 1 << " (";
                        std::cout << (seg_ori ? seg[0] : seg[1]) << " "
                            << (seg_ori ? seg[1] : seg[0]) << ") ";
                    }
                    std::cout << std::endl;
                }
            }
        }

        // Add regular attribute
        envelope.template create_attribute<uint8_t>(
            "regular",
            lagrange::AttributeElement::Facet,
            lagrange::AttributeUsage::Scalar,
            1
        );
        auto regular_values = attribute_vector_ref<uint8_t>(envelope, "regular");
        for (size_t i=0; i<num_cycles; i++) {
            regular_values[i] = isocontour.is_cycle_regular(i) ? 1 : 0;
        }
    }
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
    stopperTime = std::chrono::high_resolution_clock::now();
    auto surface_2_end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    std::cout << "Surfacing time: " << (surface_2_end - surface_1_end) * 1e-6 << " seconds (Second marching)" << std::endl;

    lagrange::io::save_mesh(output_path + "/envelope.obj", envelope);
    envelope.initialize_edges();
    //if (!lagrange::is_manifold(envelope)) {
    //    throw std::runtime_error("ERROR: lagrange envelope is not manifold");
    //}
    lagrange::triangulate_polygonal_facets(envelope);

    
#if SAVE_CONTOUR
    
    /// Mathematica isosurfacing output:
    std::vector<std::array<double, 3>> verts_math;
    std::vector<std::array<size_t, 4>> simps_math;
    std::vector<std::vector<double>> time_math;
    std::vector<std::vector<double>> values_math;
    convert_4d_grid_col(grid, vertexMap,
                        verts_math,
                        simps_math,
                        time_math,
                        values_math);
    std::string column_iso_file = "column_iso.json";
    {
        if (std::filesystem::exists(column_iso_file.c_str())){
            std::filesystem::remove(column_iso_file.c_str());
        }
        using json = nlohmann::json;
        std::ofstream fout(column_iso_file.c_str(),std::ios::app);
        json jOut;
        jOut.push_back(json(verts_math));
        jOut.push_back(json(simps_math));
        jOut.push_back(json(time_math));
        jOut.push_back(json(values_math));
        fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
        fout.close();
    }
    /// End of Mathematica output
#endif

    auto sweep_arrangement = compute_envelope_arrangement(envelope);
    auto sweep_surface = extract_sweep_surface_from_arrangement(sweep_arrangement);
    lagrange::io::save_mesh(output_path + "/sweep_surface.msh", sweep_surface);
    lagrange::io::save_mesh(output_path + "/arrangement.msh", sweep_arrangement);
    mtet::save_mesh(output_path + "/tet_grid.msh", grid);
    return 0;
}
