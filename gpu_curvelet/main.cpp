
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
#include "gpu_curvelet.hpp"
#include "timer.hpp"

static GPUPreprocessOutput parse_preprocess_layout(const std::string &mode)
{
    if (mode == "edge-look-list" || mode == "edgelooklist" || mode == "legacy") {
        return GPUPreprocessOutput::EdgeLookList;
    }
    if (mode == "both" || mode == "compare") {
        return GPUPreprocessOutput::Both;
    }
    return GPUPreprocessOutput::NeighborCSR;
}

static GPUNeighborCSRStrategy parse_csr_strategy(const std::string &mode)
{
    if (mode == "two-pass" || mode == "twopass" || mode == "legacy") {
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
    const float dx = static_cast<float>(params.dx);
    const float dt = static_cast<float>(params.dt_deg / 180.0) * static_cast<float>(M_PI);
    const float max_k = static_cast<float>(params.max_k);
    const unsigned curvelet_style = params.curvelet_style;
    const unsigned group_max_sz = params.group_max_sz;
    const unsigned out_type = params.out_type;
    unsigned max_LookEdgeNum = params.max_look_edge_num;
    const float sx = static_cast<float>(params.sx);
    const float st = static_cast<float>(params.st);
    const unsigned LOOK_EDGE_SLOTS = params.look_edge_slots;
    const GPUPreprocessOutput layout = parse_preprocess_layout(params.preprocess_layout);

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

    // ------------------------------------------------------------------
    //> GPU preprocessing
    GPUPreprocessConfig pre_cfg;
    pre_cfg.device_id = gpu_id;
    pre_cfg.num_edges = edge_num;
    pre_cfg.sz_edge_data = edge_data_sz;
    pre_cfg.look_slots = LOOK_EDGE_SLOTS;
    pre_cfg.neighbor_radius = 3;
    pre_cfg.rad = nrad;
    pre_cfg.copy_to_host = true;
    pre_cfg.output_layout = layout;
    pre_cfg.csr_strategy = parse_csr_strategy(params.csr_strategy);
    pre_cfg.max_candidates = params.max_candidates;

    GPUPreprocessResult pre_result;
    if (!gpu_preprocess_build(pre_cfg, TOED_edges.data(), pre_result)) {
        return false;
    }
    timer.lap("GPU preprocessing");

    const bool b_have_edge_look = (layout == GPUPreprocessOutput::EdgeLookList) || (layout == GPUPreprocessOutput::Both);
    const bool b_have_csr = (layout == GPUPreprocessOutput::NeighborCSR) || (layout == GPUPreprocessOutput::Both);

    if (b_have_csr) {
        max_LookEdgeNum = pre_result.csr.max_neighbor_degree;
        std::cout << "CSR layout ready: " << pre_result.csr.total_neighbor_pairs
                  << " anchor-neighbor pairs" << std::endl;
    }
    if (b_have_edge_look) {
        max_LookEdgeNum = pre_result.max_look_edge_num;
    }

    if (b_have_edge_look) {
        float *edgeLookList = pre_result.host_edgeLookList;
        const int edge_look_stride = pre_result.edge_look_stride;

        int edge_num_arg = edge_num;
        int edge_data_sz_arg = edge_data_sz;
        CurveletGPU CurveletGPU_obj(gpu_id, edge_num_arg, edge_data_sz_arg, edgeLookList, max_LookEdgeNum,
                                    edge_look_stride, dx, dt, sx, st, max_k, group_max_sz);
        timer.lap("CurveletGPU setup (alloc + copy edgeLookList)");

        CurveletGPU_obj.preprocessing();
        timer.lap("GPU curvelet init (H2D memcpy)");

        CurveletGPU_obj.build_curvelets_greedy();
        timer.lap("build_curvelets_greedy (host wall time)");
    } 
    else {
        std::cout << "Skipping legacy curvelet kernel (CSR-only preprocess mode)." << std::endl;
    }

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
