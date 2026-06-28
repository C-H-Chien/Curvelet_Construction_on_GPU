/*****************************************************************************
// file: main - GPU
// author: Chiang-Heng Chien
// date: 06/30/2022
//       An algorithm to form curvelet from input third order edge list using GPU
******************************************************************************/

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "data_io.hpp"
#include "param_settings.hpp"
#include "gpu_preprocess.hpp"
#include "gpu_common.hpp"
#include "timer.hpp"

static GPUNeighborCSRStrategy parse_csr_strategy(const std::string &mode)
{
    if (mode == "two-pass" || mode == "twopass") {
        return GPUNeighborCSRStrategy::TwoPass;
    }
    if (mode == "compare-csr" || mode == "compare") {
        return GPUNeighborCSRStrategy::CompareCSR;
    }
    return GPUNeighborCSRStrategy::SinglePass;
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

    StepTimer timer;
    timer.start();

    GPUPreprocessConfig pre_cfg;
    pre_cfg.device_id = gpu_id;
    pre_cfg.num_edges = edge_num;
    pre_cfg.sz_edge_data = edge_data_sz;
    pre_cfg.neighbor_radius = 3;
    pre_cfg.rad = nrad;
    pre_cfg.copy_to_host = true;
    pre_cfg.csr_strategy = parse_csr_strategy(params.csr_strategy);
    pre_cfg.max_candidates = params.max_candidates;

    GPUPreprocessResult pre_result;
    if (!gpu_preprocess_build(pre_cfg, TOED_edges.data(), pre_result)) {
        return false;
    }
    timer.lap("GPU CSR preprocessing");

    std::cout << "CSR layout ready: " << pre_result.csr.total_neighbor_pairs
              << " anchor-neighbor pairs, max degree "
              << pre_result.csr.max_neighbor_degree << std::endl;
    std::cout << "GPU curvelet kernel (CSR-native) not yet implemented; preprocess only."
              << std::endl;

    timer.summary();

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
