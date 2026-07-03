/*****************************************************************************
// file: main - GPU
// author: Chiang-Heng Chien
// date: 06/30/2022
//       An algorithm to form curvelet from input third order edge list using GPU
******************************************************************************/

#include <string>
#include <vector>
#define _USE_MATH_DEFINES // Must be before #include <cmath>
#include <iostream>
#include <cmath>

#include "data_io.hpp"
#include "param_settings.hpp"
#include "gpu_preprocess.hpp"
#include "gpu_curve_bundle_formation.hpp"
#include "gpu_common.hpp"
#include "timer.hpp"

static GPUNeighborCSRStrategy parse_csr_strategy(const std::string &mode)
{
    if (mode == "two-pass" || mode == "twopass") {
        return GPUNeighborCSRStrategy::TwoPass;
    }
    return GPUNeighborCSRStrategy::SinglePass;
}

static GPUFixedRowBuildStrategy parse_fixed_row_build(const std::string &mode)
{
    if (mode == "stage") {
        return GPUFixedRowBuildStrategy::Stage;
    }
    if (mode != "warp") {
        std::cerr << "Warning: unknown fixed-row build '" << mode
                  << "', using warp (expected: warp | stage)\n";
    }
    return GPUFixedRowBuildStrategy::Warp;
}

static GPUNeighborLayout parse_neighbor_layout(const std::string &mode)
{
    if (mode == "fixed-row" || mode == "fixedrow" || mode == "look-list" || mode == "looklist") {
        return GPUNeighborLayout::FixedRow;
    }
    if (mode != "csr") {
        std::cerr << "Warning: unknown neighbor layout '" << mode
                  << "', using csr (expected: csr | fixed-row)\n";
    }
    return GPUNeighborLayout::CSR;
}

bool run_curvelet_gpu(const std::string &out_chain_file, int gpu_id, const CurveletParams &params)
{
    const float nrad = static_cast<float>(params.nrad);
    const unsigned curvelet_style = params.curvelet_style;
    const unsigned out_type = params.out_type;

    const std::string &edge_file = params.edge_file;
    const int edge_data_sz = params.edge_data_sz;

    std::cout << "Using scalar type: float (GPU)" << std::endl;

    std::vector<float> TOED_edges;
    if (!read_TO_edges_from_file(edge_file, edge_data_sz, TOED_edges)) {
        return false;
    }
    int edge_num = static_cast<int>(TOED_edges.size() / edge_data_sz);

    cudaDeviceProp prop;
    cudacheck(cudaSetDevice(gpu_id));
    cudaGetDeviceProperties(&prop, gpu_id);
    printf("Device name: %s (Compute capability: %d.%d)\n", prop.name, prop.major, prop.minor);

    GPUPreprocessConfig pre_cfg;
    pre_cfg.device_id = gpu_id;
    pre_cfg.num_edges = edge_num;
    pre_cfg.sz_edge_data = edge_data_sz;
    pre_cfg.neighbor_radius = 3;
    pre_cfg.rad = nrad;
    pre_cfg.csr_strategy = parse_csr_strategy(params.csr_strategy);
    pre_cfg.neighbor_layout = parse_neighbor_layout(params.neighbor_layout);
    pre_cfg.max_candidates = params.max_candidates;
    pre_cfg.neighbor_count_threads = params.neighbor_count_threads;
    pre_cfg.neighbor_fill_threads = params.neighbor_fill_threads;
    pre_cfg.neighbor_stage_threads = params.neighbor_stage_threads;
    pre_cfg.neighbor_warps_per_block = params.neighbor_warps_per_block;
    pre_cfg.fixed_row_build = parse_fixed_row_build(params.fixed_row_build);

    CategoryProfiler profiler;
    profiler.set_title("preprocess");
    profiler.start();
    GPUPreprocessResult pre_result;
    if (!gpu_preprocess_build(pre_cfg, TOED_edges.data(), pre_result, &profiler)) {
        return false;
    }
    profiler.summary();

    if (pre_result.csr.layout == GPUNeighborLayout::FixedRow) {
        std::cout << "Fixed-row layout ready: " << pre_result.csr.total_neighbor_pairs
                  << " anchor-neighbor pairs, " << pre_result.csr.neighbor_slots_per_anchor
                  << " slots/anchor, max degree " << pre_result.csr.max_neighbor_degree
                  << std::endl;

        GPUCurveletConfig bundle_cfg;
        bundle_cfg.dx = static_cast<float>(params.dx);
        bundle_cfg.dt = static_cast<float>(params.dt_deg * M_PI / 180.0);
        bundle_cfg.sx = static_cast<float>(params.sx);
        bundle_cfg.st = static_cast<float>(params.st);
        bundle_cfg.max_k = static_cast<float>(params.max_k);
        bundle_cfg.sz_edge_data = edge_data_sz;
        bundle_cfg.bundle_warps_per_block = params.neighbor_warps_per_block;

        CategoryProfiler bundle_profiler;
        bundle_profiler.set_title("pairwise_curve_bundles");
        bundle_profiler.start();

        GPUCurveBundleStorage bundle_storage;
        GPUCurveletFormationResult bundle_result;
        if (!gpu_form_pairwise_bundles(bundle_cfg, pre_result.csr, bundle_storage, bundle_result, &bundle_profiler)) {
            gpu_curvelet_free_bundles(bundle_storage);
            gpu_preprocess_free(pre_result);
            return false;
        }
        bundle_profiler.summary();

        std::cout << "Pairwise curve bundles formed (fixed-row warp): "
                  << bundle_result.valid_pairs << " valid pairs" << std::endl;
        std::cout << "Chain growth / curvelet output not yet implemented on GPU." << std::endl;

        gpu_curvelet_free_bundles(bundle_storage);
    }
    else {
        std::cout << "CSR layout ready: " << pre_result.csr.total_neighbor_pairs
                  << " anchor-neighbor pairs, max number of neighbor edges per anchor = "
                  << pre_result.csr.max_neighbor_degree << std::endl;
        std::cout << "GPU curve bundle formation requires --neighbor-layout fixed-row." << std::endl;
    }

    gpu_preprocess_free(pre_result);

    (void)out_chain_file;
    (void)curvelet_style;
    (void)out_type;
    return true;
}

int main(int argc, char **argv)
{
    int gpu_id = 0;
    std::string out_file = "chain_gpu.txt";
    CurveletParams params;
    bool show_help = false;

    if (!parse_args(argc, argv, params, out_file, gpu_id, show_help)) {
        return show_help ? 0 : 1;
    }

    const bool ok = run_curvelet_gpu(out_file, gpu_id, params);
    return ok ? 0 : 1;
}
