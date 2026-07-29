#ifndef GPU_CURVELET_HPP
#define GPU_CURVELET_HPP

#include "gpu_preprocess.hpp"
#include "param_settings.hpp"

#include <cmath>

class CategoryProfiler;

struct GPUCurveBundleStorage {
    int curves_num_in_bundle_pixel = 0;
    int curves_num_in_bundle_theta = 0;
    int bundle_cells = 0;
    int slots_per_anchor = 0;
    int num_edges = 0;

    //> Slot-major within each anchor: [anchor][slot][cell]
    float *dev_bundle_min_ks = nullptr;
    float *dev_bundle_max_ks = nullptr;
    unsigned char *dev_hyp_look_edge = nullptr;
};

//> Calculate the perturbation region size of edge location and orientation
inline int get_gpu_curve_bundle_grid_size(float tolerance, float sample_step)
{
    return 2 * static_cast<int>(floorf(tolerance / sample_step + 0.5f)) + 1;
}

bool gpu_allocate_curve_bundles(
    const GPUNeighborGraph &graph,
    GPUCurveBundleStorage &storage,
    unsigned **dev_valid_pair_count);


//> Pairwise curve bundle formation on a fixed-row neighbor graph (one warp per anchor)
bool gpu_form_pairwise_bundles_main(
    const CurveletParams &params,
    const GPUNeighborGraph &graph,
    GPUCurveBundleStorage &storage,
    unsigned &valid_pairs,
    CategoryProfiler *profiler = nullptr);

void gpu_curvelet_free_bundles(GPUCurveBundleStorage &storage);

#endif // GPU_CURVELET_HPP
