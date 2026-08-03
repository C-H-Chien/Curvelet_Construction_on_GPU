#ifndef GPU_EDGE_CHAIN_GROWTH_HPP
#define GPU_EDGE_CHAIN_GROWTH_HPP

#include "gpu_curve_bundle_formation.hpp"

#include <string>
#include <vector>

class CategoryProfiler;

struct GPUCurveletChainStorage {
    int num_edges = 0;
    int slots_per_anchor = 0;
    int group_max_sz = 0;
    int chain_width = 0;       // group_max_sz + 1
    int max_per_anchor = 0;
    int max_curvelets = 0;
    int bundle_cells = 0;

    //> Warp growth shared-memory plan (set at allocate / launch time)
    int warp_smem_mode = 0;          // 0=none, 1=lane workspace in shared
    int warp_warps_per_block = 0;    // equals request unless --chain-smem-mode auto reduces it
    size_t warp_smem_bytes = 0;

    unsigned *dev_edge_chain_final = nullptr;   // max_curvelets * chain_width
    unsigned *dev_anchor_chain_count = nullptr; // num_edges
    float *dev_curvelet_info = nullptr;          // max_curvelets * CURVELET_INFO_WIDTH (row-major)

    //> Per-anchor scratch for compare / intersect / on-the-fly / dedup buffers
    float *dev_scratch_f = nullptr;
    unsigned *dev_scratch_u = nullptr;
    size_t scratch_floats_per_anchor = 0;
    size_t scratch_uints_per_anchor = 0;
};

struct GPUCurveletChainResult {
    unsigned num_curvelets = 0;
};

inline constexpr int GPU_CURVELET_INFO_WIDTH = 10;

bool gpu_allocate_edge_chains(
    const GPUNeighborGraph &graph,
    const GPUCurveBundleStorage &bundles,
    int group_max_sz,
    int chain_warps_per_block,
    const std::string &chain_smem_mode,
    GPUCurveletChainStorage &storage);

void gpu_curvelet_free_chains(GPUCurveletChainStorage &storage);

//> Grow edge chains by cumulative curve-bundle intersection (style-2: anchor-leading, both directions).
//> Direction filtering is applied here against pairwise hyp flags from bundle formation.
//> Uses warp-per-anchor growth (phase 1) then thread-per-anchor dedup/record (phase 2).
bool gpu_grow_edge_chains_main(
    const CurveletParams &params,
    const GPUNeighborGraph &graph,
    const GPUCurveBundleStorage &bundles,
    GPUCurveletChainStorage &storage,
    GPUCurveletChainResult &result,
    CategoryProfiler *profiler = nullptr);

//> Compact per-anchor chain rows into a contiguous host buffer (1-based edge IDs, CPU-compatible).
bool gpu_download_compact_chains(
    const GPUCurveletChainStorage &storage,
    std::vector<int> &host_chains,
    std::vector<float> &host_info,
    unsigned &num_curvelets);

#endif // GPU_EDGE_CHAIN_GROWTH_HPP
