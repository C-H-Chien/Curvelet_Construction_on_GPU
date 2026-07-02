#ifndef GPU_CURVELET_HPP
#define GPU_CURVELET_HPP

#include "gpu_preprocess.hpp"

#include <cmath>

class CategoryProfiler;

struct GPUCurveletConfig {
    float dx = 0.4f;
    float dt = 0.261799f;  // 15 deg
    float sx = 0.1f;
    float st = 0.08f;
    float max_k = 0.3f;
    int sz_edge_data = 4;
    int bundle_warps_per_block = 4;
};

struct GPUCurveBundleStorage {
    int curves_num_in_bundle_pixel = 0;
    int curves_num_in_bundle_theta = 0;
    int bundle_cells = 0;
    int slots_per_anchor = 0;
    int num_edges = 0;

    float *dev_bundle_min_ks = nullptr;
    float *dev_bundle_max_ks = nullptr;
    unsigned char *dev_hyp_look_edge = nullptr;
};

struct GPUCurveletFormationResult {
    unsigned valid_pairs = 0;
};

//> Calculate the perturbation region size of edge location and orientation
inline int gpu_curve_bundle_grid_size(float tolerance, float sample_step)
{
    return 2 * static_cast<int>(floorf(tolerance / sample_step + 0.5f)) + 1;
}

bool gpu_allocate_curve_bundles(
    const GPUNeighborGraph &graph,
    GPUCurveBundleStorage &storage,
    unsigned **dev_valid_pair_count);
void gpu_curvelet_free_bundles(GPUCurveBundleStorage &storage);

//> Pairwise curve bundle formation on a fixed-row neighbor graph (one warp per anchor).
bool gpu_form_pairwise_bundles(
    const GPUCurveletConfig &cfg,
    const GPUNeighborGraph &graph,
    GPUCurveBundleStorage &storage,
    GPUCurveletFormationResult &result,
    CategoryProfiler *profiler = nullptr);

#endif // GPU_CURVELET_HPP
