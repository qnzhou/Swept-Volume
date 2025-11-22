//
//  post_processing.h
//  general_sweep
//
//  Created by Yiwen Ju on 5/13/25.
//

#ifndef post_processing_h
#define post_processing_h
#include <numeric>   // for std::iota
#include <utility>   // for std::pair
#include <limits>    // for std::numeric_limits
#include <chrono>
#include <iostream>
#include <algorithm>

#include <ankerl/unordered_dense.h>

#include <arrangement/Arrangement.h>

#include <lagrange/SurfaceMesh.h>
#include <lagrange/combine_meshes.h>
#include <lagrange/compute_area.h>
#include <lagrange/separate_by_components.h>
#include <lagrange/views.h>

const double threshold = 1e-5;
const int faceThreshold = 200;

struct PairHasher {
    std::size_t operator()(const std::pair<int,int>& p) const noexcept {
        auto a = static_cast<uint32_t>(p.first);
        auto b = static_cast<uint32_t>(p.second);
        if (a > b) std::swap(a, b);
        uint64_t high = static_cast<uint64_t>(a) << 32;
        uint64_t low  = static_cast<uint64_t>(b);
        uint64_t packed = high | low;
        return static_cast<std::size_t>(packed);
    }
};
struct PairEqual {
    bool operator()(std::pair<int,int> x,
                    std::pair<int,int> y) const noexcept
    {
        if (x.first > x.second) std::swap(x.first,x.second);
        if (y.first > y.second) std::swap(y.first,y.second);
        return x.first==y.first && x.second==y.second;
    }
};
constexpr std::pair<int,int> INVALID_PATCH = { -1, -1 };

std::pair<std::vector<bool>, std::vector<int>>
computeValidPatch(int cellNum,
                  std::vector<double>& volInfo,
                  std::vector<int>& faceCountInfo,
                  const arrangement::MatrixIr& cellData)
{
    std::vector<std::vector<bool>> adj(cellNum, std::vector<bool>(cellNum, false));
    std::vector<bool> valid(cellNum, false), merged(cellNum, false);
    for(int i = 0; i < cellNum; ++i) {
        if (volInfo[i] > threshold && faceCountInfo[i] > faceThreshold) valid[i] = true;
    }
    // build adjacency from cellData
    for( size_t e = 0; e < cellData.rows(); e++) {
        int u = cellData(e, 0), v = cellData(e, 1);
            adj[u][v] = adj[v][u] = true;
    }
    std::vector<int> owner(cellNum);
    for(int i = 0; i < cellNum; ++i) owner[i] = i;
    int cluster_size = cellNum;
    auto terminate = [&](){
        for(int i = 0; i < cluster_size ; i++)
            if (!(valid[i] || merged[i])) return false;
        return true;
    };

    while(!terminate()) {
        int mergeCell = 0;
        while(mergeCell < cluster_size && (valid[mergeCell] || merged[mergeCell]))
            ++mergeCell;
            // Do not need to consider cells that are either valid or merged.

        int target = -1;
        double target_vol = std::numeric_limits<double>::lowest();
        for (int clusterIt = 0; clusterIt < cluster_size; ++clusterIt) {
            if (!merged[clusterIt] && adj[mergeCell][clusterIt] && volInfo[clusterIt] > target_vol){
                target_vol = volInfo[clusterIt];
                target = clusterIt;
            }
        }
        assert(target < cluster_size);
        if (target >= cluster_size) break;    // nothing left to merge
        // merging two clusters
        cluster_size ++;
        // update ownership
        for(int j = 0; j < cellNum; ++j){
          if(owner[j] == target || owner[j] == mergeCell)
            owner[j] = cluster_size - 1;
        }
        merged[target] = merged[mergeCell] = true;
        merged.push_back(false);
        valid.push_back(valid[target] || valid[mergeCell]);
        volInfo.push_back(std::max(volInfo[target], volInfo[mergeCell]));
        std::vector<bool> new_cluster(cluster_size, false);
        for(int i = 0; i < cluster_size - 1; ++i) {
            new_cluster[i] = adj[target][i] || adj[mergeCell][i];
        }
        new_cluster[target] = new_cluster[mergeCell] = false;
        for (int i = 0; i < cluster_size - 1; ++i){
            adj[i].push_back(new_cluster[i]);
        }
        adj.push_back(new_cluster);
    }
    std::vector<bool> validPatch(cellData.rows(), false);
    for(int i = 0; i < cellData.rows(); ++i) {
        int u = cellData(i, 0),
            v = cellData(i, 1);
        if (owner[u] != owner[v]) {
            validPatch[i] = true;
        }
    }
    return {validPatch, owner};
}

template <typename Scalar, typename Index>
lagrange::SurfaceMesh<Scalar, Index> compute_swept_volume_from_envelope(
        const lagrange::SurfaceMesh<Scalar, Index>& envelope) {

    using Point = Eigen::Matrix<Scalar, 3, 1>;

    // Compute arrangement
    auto V = vertex_view(envelope).template cast<double>();
    auto F = facet_view(envelope).template cast<int>();
    auto T = attribute_vector_view<Scalar>(envelope, "time");
    arrangement::VectorI face_labels = Eigen::VectorXi::LinSpaced(F.rows(), 0, F.rows() - 1);
    auto engine = arrangement::Arrangement::create_mesh_arrangement(V, F, face_labels);
    engine->run();

    lagrange::SurfaceMesh<Scalar, Index> output_mesh;

    const auto& cell_data       = engine->get_cells();           // (#facets x 2) array
    const auto& patches         = engine->get_patches();         // list of facet indices
    const auto& parent_facets   = engine->get_out_face_labels();   // size = #facets
    const auto& winding_number  = engine->get_winding_number();  // (#facets x 2) array

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
    for (size_t fid=0; fid < num_facets; fid++) {
        Index c0 = static_cast<Index>(cell_data(patches[fid], 0)); // Cell on the positive side
        Index c1 = static_cast<Index>(cell_data(patches[fid], 1)); // Cell on the negative side

        cell_graph[c0].insert(c1);
        cell_graph[c1].insert(c0);

        Eigen::Matrix<int, 1, 3> f = arrangement_faces.row(fid);
        Point p0 = out_vertices.row(f[0]).template cast<Scalar>();
        Point p1 = out_vertices.row(f[1]).template cast<Scalar>();
        Point p2 = out_vertices.row(f[2]).template cast<Scalar>();
        Scalar vol = p0.dot( (p1).cross(p2) );

        cell_volumes[c0] += vol / 6;
        cell_volumes[c1] -= vol / 6;
    }
    std::vector<bool> cell_is_isolated(num_cells, false);
    for (size_t cid=0; cid < num_cells; cid++) {
        if (cell_graph[cid].size() <= 1) {
            // Cell is completely embedded in another cell
            // It can be safely removed without affecting other cells and their winding numbers
            cell_is_isolated[cid] = true;
        }
    }

    output_mesh.add_vertices(out_vertices.rows(),
            {out_vertices.data(), static_cast<size_t>(out_vertices.size())});

    size_t num_valid_facets = 0;
    std::vector<int8_t> valid(num_facets, 0);
    for (size_t fid=0; fid < num_facets; fid++) {
        int w0 = winding_number(fid, 0); // Winding number on the positive side
        int w1 = winding_number(fid, 1); // Winding number on the negative side
        int c0 = cell_data(patches[fid], 0); // Cell on the positive side
        int c1 = cell_data(patches[fid], 1); // Cell on the negative side

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
    for (size_t fid=0; fid < num_facets; fid++) {
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
        "time",
        lagrange::AttributeElement::Corner,
        lagrange::AttributeUsage::Scalar,
        1
    );
    auto time_values = attribute_vector_ref<Scalar>(output_mesh, "time");
    count = 0;
    // TODO: parallalize this loop
    for (size_t fid=0; fid < num_facets; fid++) {
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

template <typename Scalar, typename Index>
lagrange::SurfaceMesh<Scalar, Index> filter_tiny_components(
    lagrange::SurfaceMesh<Scalar, Index>& mesh,
    Scalar time_threshold)
{
    lagrange::SeparateByComponentsOptions options;
    options.map_attributes = true;
    auto comps = lagrange::separate_by_components(mesh, options);

    auto compute_volume = [](const lagrange::SurfaceMesh<Scalar, Index>& comp) {
        Scalar comp_vol = 0;
        auto V = vertex_view(comp);
        auto F = facet_view(comp);
        auto num_facets = comp.get_num_facets();
        for (Index fid=0; fid < num_facets; fid++) {
            Eigen::Matrix<Scalar, 1, 3> p0 = V.row(F(fid, 0));
            Eigen::Matrix<Scalar, 1, 3> p1 = V.row(F(fid, 1));
            Eigen::Matrix<Scalar, 1, 3> p2 = V.row(F(fid, 2));
            Scalar vol = p0.dot( (p1).cross(p2) );
            comp_vol += vol;
        }
        return comp_vol / 6;
    };

    auto compute_time_span = [](const lagrange::SurfaceMesh<Scalar, Index>& comp) {
        auto time_values = attribute_vector_view<Scalar>(comp, "time");
        return time_values.maxCoeff() - time_values.minCoeff();
    };

    std::erase_if(comps, [&](const auto& comp){
        Scalar vol = compute_volume(comp);
        Scalar area = lagrange::compute_mesh_area(comp);
        Scalar time_span = compute_time_span(comp);
        std::cout << "Volume: " << vol << "  Area: " << area << " time span: " << time_span << std::endl;
        return time_span < time_threshold;
    });

    return lagrange::combine_meshes<Scalar, Index>(
            std::span<lagrange::SurfaceMesh<Scalar, Index>>(comps.data(), comps.size()));
}


/**
 * Compute the surface of the sweep volume.
 *
 * @param vertices Input vertices of the sweep generator.
 * @param faces Input faces of the sweep generator.
 * @param out_vertices Output vertices of the sweep volume surface.
 * @param out_faces Output faces of the sweep volume surface.
 */
void compute_sweep_volume(const arrangement::MatrixFr& vertices, const arrangement::MatrixIr& faces,
                          arrangement::MatrixFr& out_vertices, std::vector<arrangement::MatrixIr>& out_faces, std::string output_path) {
    // assume the input faces have already been oriented based on second order time derivatives
    // then the sweep volume consists of arrangement cells with positive winding number
    
    using Clock = std::chrono::high_resolution_clock;
    using Milliseconds = std::chrono::milliseconds;

    Clock::time_point t0 = Clock::now();

    // Initialize face labels
    arrangement::VectorI face_labels = Eigen::VectorXi::LinSpaced(faces.rows(), 0, faces.rows() - 1);
    // create a mesh arrangement engine.
    auto engine = arrangement::Arrangement::create_mesh_arrangement(vertices, faces, face_labels);
    // Alternatively, create a fast mesh arrangement engine.
    // face_labels.setZero();  // fast arrangement cannot use too large face labels
    // auto engine = arrangement::Arrangement::create_fast_arrangement(vertices, faces, face_labels);
    // or, create a Geogram arrangement engine
    // the cell ids associated with patches can often be MAX_INT64.
    // auto engine = arrangement::Arrangement::create_geogram_arrangement(vertices, faces, face_labels);
    
// The following is based on James' arrangement code in Python: https://github.com/qnzhou/arrangement-benchmark/blob/main/python/arrangement/__main__.py
    engine->run();

    Clock::time_point t1 = Clock::now();
    
    const auto& cell_data       = engine->get_cells();           // (#facets x 2) array
    const auto& patches         = engine->get_patches();         // list of facet indices
    const auto& parent_facets   = engine->get_out_face_labels();   // size = #facets
    const auto& winding_number  = engine->get_winding_number();  // (#facets x 2) array
    out_vertices = engine->get_vertices();
    const auto& arrangement_faces = engine->get_faces();
    int num_cells = engine->get_num_cells();
    int num_patches = engine->get_num_patches();
    int num_facets = arrangement_faces.rows();
    assert(patches.size() == num_facets);
    std::vector<int> wind_list(num_cells, std::numeric_limits<int>::max());
    std::vector<double> volInfo(num_cells, 0);
    std::vector<int> faceCountInfo(num_cells, 0);
    std::vector<size_t> cellIt;

    for (size_t facet_idx = 0; facet_idx < patches.size(); facet_idx++) {
        int patch_idx = patches[facet_idx];
        int c0 = cell_data(patch_idx, 0); // Cell on the positive side
        int c1 = cell_data(patch_idx, 1); // Cell on the negative side
        int j = arrangement_faces(facet_idx, 0), k = arrangement_faces(facet_idx, 1), l = arrangement_faces(facet_idx, 2);
        Eigen::Vector3d vj = out_vertices.row(j);
        Eigen::Vector3d vk = out_vertices.row(k);
        Eigen::Vector3d vl = out_vertices.row(l);
        double vol = vj.dot( vk.cross(vl) );
        volInfo[c0] -= vol / 6;
        volInfo[c1] += vol / 6;
        faceCountInfo[c0] += 1;
        faceCountInfo[c1] += 1;
        if (wind_list[c0] == std::numeric_limits<int>::max()) {
            wind_list[c0] = winding_number(facet_idx, 0);
        } else {
            assert(wind_list[c0] == winding_number(facet_idx, 0));
        }
        if (wind_list[c1] == std::numeric_limits<int>::max()) {
            wind_list[c1] = winding_number(facet_idx, 1);
        } else {
            assert(wind_list[c1] == winding_number(facet_idx, 1));
        }
    }

    std::vector<bool> valid(num_cells, false);
    int valid_num = 0;
    std::cout << "valid 0-winding cell iter: ";
    for (size_t i = 0; i < num_cells; ++i) {
        volInfo[i] = std::abs(volInfo[i]) / 6.0;
        if (volInfo[i] > threshold && faceCountInfo[i] > faceThreshold) {
            valid[i] = true;
            valid_num++;
        }
        if (valid[i] && wind_list[i] == 0){
            cellIt.push_back(i);
            std::cout << i << " ";
        }
    }
    std::cout << std::endl;
    std::cout << "number of valid cells: " << valid_num << std::endl;
    
    // Pruning of error cells
    const auto& prunedInfo2 = computeValidPatch(num_cells, volInfo, faceCountInfo, cell_data);
    std::vector<bool> valid_patchInd = prunedInfo2.first;
    std::vector<int> cell_merge_map = prunedInfo2.second;
//    std::vector<bool> valid_patchInd (cell_data.size(), true);
//    std::vector<int> cell_merge_map (num_cells, -1);
//    for (int i = 0; i < num_cells; i++){
//        cell_merge_map[i] = i;
//    }
    
    std::vector<std::array<std::array<double, 3>, 2>> feature_lines;
    std::vector<std::array<double, 3>> corners;
    std::unordered_map<int, std::pair<int, int>> patchIdMap;
    for (int i = 0; i < valid_patchInd.size(); i++){
        if (valid_patchInd[i]){
            int c0 = cell_data(i, 0);
            int c1 = cell_data(i, 1);
            patchIdMap[i].first = cell_merge_map[c0];
            patchIdMap[i].second = cell_merge_map[c1];
        }
    }
    
    // New detection method to check if multiple patches meet (valences information is also needed after testing)
    std::unordered_map<std::pair<int, int>, std::pair<int, int>, PairHasher, PairEqual> edge_patch_map;
    std::unordered_map<std::pair<int, int>, int, PairHasher, PairEqual> edge_valence_map;
    std::unordered_map<int, int> vert_valence;
    auto push_vert_valance = [&](int v){
        vert_valence[v] ++;
    };
    auto push_edge_patch = [&](std::pair<int, int> edge, int patchId){
        int c0 = cell_data(patchId, 0), c1 = cell_data(patchId, 1);
        if (c0 > c1){
            std::swap(c0, c1);
        }
        edge_patch_map[edge] = std::pair<int, int>{c0, c1};
        edge_valence_map[edge]++;
    };
    auto check_edge_patch = [&](std::pair<int, int> edge, int patchId){
        edge_valence_map[edge]++;
        auto prev_patch = edge_patch_map[edge];
        if (prev_patch != INVALID_PATCH){
            PairHasher hasher;
            if (hasher(prev_patch) != hasher({cell_data(patchId, 0), cell_data(patchId, 1)})){
                edge_patch_map[edge] = INVALID_PATCH;
            }
        }
    };
    auto parse_edge_patch = [&] (std::pair<int, int> edge){
        if (edge_patch_map[edge] == INVALID_PATCH && edge_valence_map[edge] > 2){
            edge_valence_map[edge] = 2;
            push_vert_valance(edge.first);
            push_vert_valance(edge.second);
            feature_lines.push_back(std::array<std::array<double, 3>, 2>{std::array<double, 3>{out_vertices(edge.first, 0), out_vertices(edge.first, 1), out_vertices(edge.first, 2)}, std::array<double, 3>{out_vertices(edge.second, 0), out_vertices(edge.second, 1), out_vertices(edge.second, 2)}});
        }
    };

    for (size_t i = 0; i < arrangement_faces.rows(); i++){
        auto& p = patches(i);
        if (valid_patchInd[p]){
            const auto& face = arrangement_faces.row(i);
            std::array<std::pair<int, int>, 3> edge =
                {std::pair<int, int>{face[0], face[1]},
                std::pair<int, int>{face[0], face[2]},
                std::pair<int, int>{face[1], face[2]}};
            for (auto& e: edge){
                if (edge_patch_map.contains(e)){
                    check_edge_patch(e, p);
                }else{
                    push_edge_patch(e, p);
                }
            }
        }
    }
    
    std::vector<std::pair<size_t, bool>> patch_list;
    for (size_t i = 0; i < cellIt.size(); ++i){
        auto current_cell = cell_merge_map[cellIt[i]];
//        auto current_cell = cellIt[i];
        arrangement::MatrixIr patch_faces;
        patch_faces.resize(arrangement_faces.rows(), 3);
        size_t face_count = 0;
        for (size_t f = 0; f < num_facets; ++f){
            auto patch_id = patches(f);
            auto c0 = cell_merge_map[cell_data(patch_id, 0)], c1 = cell_merge_map[cell_data(patch_id, 1)];
//            auto c0 = cell_data(patch_id, 0), c1 = cell_data(patch_id, 1);
            if ((c0 == current_cell && c1 != current_cell) || (c1 == current_cell && c0 != current_cell)){
                if (c0 == current_cell){
                    patch_faces.row(face_count) = (arrangement_faces.row(f));
                }else{
                    patch_faces.row(face_count) = (arrangement_faces.row(f).reverse());
                }
                const auto& face = arrangement_faces.row(f);
                std::array<std::pair<int, int>, 3> edge =
                    {std::pair<int, int>{face[0], face[1]},
                    std::pair<int, int>{face[0], face[2]},
                    std::pair<int, int>{face[1], face[2]}};
                for (auto& e: edge){
                    parse_edge_patch(e);
                }
                face_count++;
            }
        }
        patch_faces.conservativeResize(face_count, 3);
        out_faces.emplace_back(patch_faces);
    }
    int count = 0;
    for (int i = 0; i < out_vertices.rows(); i++){
        if (vert_valence[i] != 2 && vert_valence[i] > 0){
            corners.push_back({out_vertices(i, 0), out_vertices(i, 1), out_vertices(i, 2)});
        }
    }

    Clock::time_point t2 = Clock::now();

    std::cout << "# Verts: " << engine->get_vertices().rows() << std::endl;
    std::cout << "# Tris: " << arrangement_faces.rows() << std::endl;
    std::cout << "patches: " << num_patches << std::endl;
    std::cout << "cells: " << num_cells << std::endl;
    std::string feature_lines_file = output_path + "/features.json";
    if (std::filesystem::exists(feature_lines_file.c_str())){
        std::filesystem::remove(feature_lines_file.c_str());
    }
    using json = nlohmann::json;
    std::ofstream fout(feature_lines_file.c_str(),std::ios::app);
    json jOut;
    jOut["chains"] = feature_lines;
    jOut["joints"] = corners;
    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
    fout.close();

    Clock::time_point t3 = Clock::now();

    auto arrangement_time = std::chrono::duration_cast<Milliseconds>(t1 - t0).count();
    auto post_process_time = std::chrono::duration_cast<Milliseconds>(t2 - t1).count();
    auto feature_output_time = std::chrono::duration_cast<Milliseconds>(t3 - t2).count();
    std::cout << "Arrangement time: " << arrangement_time * 1e-3 << " seconds" << std::endl;
    std::cout << "Post-process time: " << post_process_time * 1e-3 << " seconds" << std::endl;
    std::cout << "Feature output time: " << feature_output_time * 1e-3 << " seconds" << std::endl;
}
#endif /* post_processing_h */
