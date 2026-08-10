/*****************************************************************************
// file: main - GPU
// author: Chiang-Heng Chien
// date: 06/30/2022
//       An algorithm to form curvelet from input third order edge list using GPU
******************************************************************************/

#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "data_io.hpp"
#include "param_settings.hpp"
#include "gpu_preprocess.hpp"
#include "gpu_curve_bundle_formation.hpp"
#include "gpu_edge_chain_growth.hpp"
#include "gpu_common.hpp"
#include "timer.hpp"

bool run_curvelet_gpu(const std::string &out_chain_file, int gpu_id, CurveletParams &params)
{
    const unsigned curvelet_style = params.curvelet_style;
    const unsigned out_type = params.out_type;

    const std::string &edge_file = params.edge_file;
    const int edge_data_sz = params.edge_data_sz;

    if (params.chain_smem_mode != "auto" && params.chain_smem_mode != "none" &&
        params.chain_smem_mode != "lane") {
        std::cerr << "Warning: unknown --chain-smem-mode '" << params.chain_smem_mode << "', using auto (expected: auto/none/lane)\n";
        params.chain_smem_mode = "auto";
    }

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

    CategoryProfiler profiler;
    profiler.set_title("preprocess");
    profiler.start();
    GPUNeighborGraph neighbor_graph;
    if (!gpu_preprocess_build(params, gpu_id, edge_num, TOED_edges.data(), neighbor_graph, &profiler)) {
        return false;
    }
    profiler.summary();

    if (neighbor_graph.layout == "fixed-row") {
        std::cout << "Fixed-row layout ready: " << neighbor_graph.total_neighbor_pairs
                  << " anchor-neighbor pairs, " << neighbor_graph.neighbor_slots_per_anchor
                  << " slots/anchor, max number of neighbors = " << neighbor_graph.max_num_of_neighbors
                  << std::endl;

        CategoryProfiler bundle_profiler;
        bundle_profiler.set_title("pairwise_curve_bundles");
        bundle_profiler.start();

        //> ================== Pairwise Curve Bundle Formation ==================
        GPUCurveBundleStorage bundle_storage;
        unsigned valid_pairs = 0;
        if (!gpu_form_pairwise_bundles_main(params, neighbor_graph, bundle_storage, valid_pairs, &bundle_profiler)) {
            gpu_curvelet_free_bundles(bundle_storage);
            gpu_preprocess_free(neighbor_graph);
            return false;
        }
        bundle_profiler.summary();
        std::cout << "Pairwise curve bundles formed (fixed-row warp): " << valid_pairs << " valid pairs" << std::endl;

        //> ================== Edge Chain Growth ==================
        CategoryProfiler chain_profiler;
        chain_profiler.set_title("grow_edge_chains");
        chain_profiler.start();

        GPUCurveletChainStorage chain_storage;
        unsigned num_curvelets = 0;
        if (!gpu_grow_edge_chains_main(params, neighbor_graph, bundle_storage, chain_storage, num_curvelets, &chain_profiler)) {
            gpu_curvelet_free_chains(chain_storage);
            gpu_curvelet_free_bundles(bundle_storage);
            gpu_preprocess_free(neighbor_graph);
            return false;
        }
        chain_profiler.summary();
        std::cout << "Edge chains grown: " << num_curvelets << " curvelets" << std::endl;

        std::vector<int> host_chains;
        std::vector<float> host_info;
        if (!gpu_download_compact_chains(chain_storage, host_chains, host_info, num_curvelets)) {
            gpu_curvelet_free_chains(chain_storage);
            gpu_curvelet_free_bundles(bundle_storage);
            gpu_preprocess_free(neighbor_graph);
            return false;
        }

        const unsigned out_w = static_cast<unsigned>(chain_storage.chain_width);
        std::cout << "(out_h, out_w) = (" << num_curvelets << ", " << out_w << ")" << std::endl;
        write_int_array_to_file(out_chain_file, host_chains.data(),
                                static_cast<int>(num_curvelets), static_cast<int>(out_w));

        std::vector<double> host_info_d(host_info.begin(), host_info.end());
        write_double_array_to_file(chain_to_info_filename(out_chain_file), host_info_d.data(),
                                   static_cast<int>(num_curvelets), GPU_CURVELET_INFO_WIDTH);

        gpu_curvelet_free_chains(chain_storage);
        gpu_curvelet_free_bundles(bundle_storage);
    }
    else {
        std::cout << "CSR layout ready: " << neighbor_graph.total_neighbor_pairs
                  << " anchor-neighbor pairs, max number of neighbor edges per anchor = "
                  << neighbor_graph.max_num_of_neighbors << std::endl;
        std::cout << "GPU curve bundle formation requires --neighbor-layout fixed-row." << std::endl;
    }

    gpu_preprocess_free(neighbor_graph);

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
