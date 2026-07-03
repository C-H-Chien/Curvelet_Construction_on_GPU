#include "gpu_curve_bundle_formation.hpp"
#include "gpu_common.hpp"
#include "timer.hpp"
#include "device/gpu_curve_bundle_formation_kernels.cuh"

#include <cmath>
#include <cstdio>

namespace {

void profile_lap(CategoryProfiler *profiler, TimerCategory cat, const char *detail)
{
    if (profiler != nullptr) {
        profiler->lap(cat, detail);
    }
}

void launch_curve_bundle_formation_kernel(
    const GPUCurveletConfig &cfg,
    const GPUNeighborGraph &graph,
    const GPUCurveBundleStorage &storage,
    unsigned *dev_valid_pair_count)
{
    //> Sanity check
    const int warps_per_block = cfg.bundle_warps_per_block;
    if (warps_per_block <= 0 || warps_per_block > 32) {
        fprintf(stderr, "launch_curve_bundle_formation_kernel: invalid bundle_warps_per_block=%d (use 1..32)\n", warps_per_block);
        return;
    }

    const int threads_per_block = warps_per_block * 32;
    const int num_blocks = div_up(graph.num_edges, warps_per_block);
    const size_t shmem_bytes = static_cast<size_t>(warps_per_block) * 2 * sizeof(float);

    curve_bundle_formation_kernel<<<num_blocks, threads_per_block, shmem_bytes>>>(
        graph.num_edges,
        storage.slots_per_anchor,
        storage.bundle_cells,
        storage.curves_num_in_bundle_pixel,
        storage.curves_num_in_bundle_theta,
        warps_per_block,
        cfg.dx, cfg.dt, cfg.sx, cfg.st, cfg.max_k,
        cfg.sz_edge_data,
        graph.dev_edges,
        graph.dev_neighbor_list,
        graph.dev_neighbor_counts,
        storage.dev_bundle_min_ks,
        storage.dev_bundle_max_ks,
        storage.dev_hyp_look_edge,
        dev_valid_pair_count);
    cudacheck(cudaGetLastError());
}

} // namespace

bool gpu_allocate_curve_bundles(
    const GPUNeighborGraph &graph,
    GPUCurveBundleStorage &storage,
    unsigned **dev_valid_pair_count)
{
    //> Sanity check
    if (graph.dev_neighbor_list == nullptr || graph.dev_neighbor_counts == nullptr) {
        fprintf(stderr, "gpu_allocate_curve_bundles: missing fixed-row neighbor buffers\n");
        return false;
    }
    if (dev_valid_pair_count == nullptr) {
        fprintf(stderr, "gpu_allocate_curve_bundles: null valid-pair counter output pointer\n");
        return false;
    }

    //> Records layout metadata
    storage.num_edges        = graph.num_edges;
    storage.slots_per_anchor = graph.neighbor_slots_per_anchor;
    storage.bundle_cells     = storage.curves_num_in_bundle_pixel * storage.curves_num_in_bundle_theta;

    if (storage.dev_bundle_min_ks == nullptr) {
        //> GPU memory size:
        //  (i)  bundle_words is the min/max curvature bounds per bundle cell
        //  (ii) hyp_words is the hypothesis lookup table for each anchor-neighbor pair (1 byte per curve bundle hypothesis for each anchor-neighbor pair)
        //       This can be 1-bit per curve nbindle but a plain unsigned char array gives coalesced 1-byte stores per thread.
        const size_t bundle_words = static_cast<size_t>(graph.num_edges) * static_cast<size_t>(storage.slots_per_anchor) * static_cast<size_t>(storage.bundle_cells);
        const size_t hyp_words    = static_cast<size_t>(graph.num_edges) * static_cast<size_t>(storage.slots_per_anchor);

        cudacheck(cudaMalloc(&storage.dev_bundle_min_ks, bundle_words * sizeof(float)));
        cudacheck(cudaMalloc(&storage.dev_bundle_max_ks, bundle_words * sizeof(float)));
        cudacheck(cudaMalloc(&storage.dev_hyp_look_edge, hyp_words    * sizeof(unsigned char)));
    }

    cudacheck(cudaMalloc(dev_valid_pair_count, sizeof(unsigned)));
    cudacheck(cudaMemset(*dev_valid_pair_count, 0, sizeof(unsigned)));
    return true;
}

void gpu_curvelet_free_bundles(GPUCurveBundleStorage &storage)
{
    if (storage.dev_bundle_min_ks != nullptr) {
        cudaFree(storage.dev_bundle_min_ks);
        storage.dev_bundle_min_ks = nullptr;
    }
    if (storage.dev_bundle_max_ks != nullptr) {
        cudaFree(storage.dev_bundle_max_ks);
        storage.dev_bundle_max_ks = nullptr;
    }
    if (storage.dev_hyp_look_edge != nullptr) {
        cudaFree(storage.dev_hyp_look_edge);
        storage.dev_hyp_look_edge = nullptr;
    }
}

bool gpu_form_pairwise_bundles(
    const GPUCurveletConfig &cfg,
    const GPUNeighborGraph &graph,
    GPUCurveBundleStorage &storage,
    GPUCurveletFormationResult &result,
    CategoryProfiler *profiler)
{
    if (graph.layout != GPUNeighborLayout::FixedRow) {
        fprintf(stderr, "gpu_form_pairwise_bundles: requires fixed-row neighbor layout\n");
        return false;
    }

    storage.curves_num_in_bundle_pixel = gpu_curve_bundle_grid_size(cfg.dx, cfg.sx);
    storage.curves_num_in_bundle_theta = gpu_curve_bundle_grid_size(cfg.dt, cfg.st);
    storage.slots_per_anchor           = graph.neighbor_slots_per_anchor;
    storage.num_edges                  = graph.num_edges;
    storage.bundle_cells               = storage.curves_num_in_bundle_pixel * storage.curves_num_in_bundle_theta;

    unsigned *dev_valid_pair_count = nullptr;

    //> Allocate the GPU memory
    const bool first_bundle_alloc = (storage.dev_bundle_min_ks == nullptr);
    if (!gpu_allocate_curve_bundles(graph, storage, &dev_valid_pair_count)) {
        return false;
    }
    profile_lap(profiler, TimerCategory::MemoryAlloc, first_bundle_alloc ? "curve bundle buffers and valid-pair counter" : "bundle valid-pair counter");

    //> Launch the curve bundle formation kernel
    launch_curve_bundle_formation_kernel(cfg, graph, storage, dev_valid_pair_count);
    profile_lap(profiler, TimerCategory::Kernel, "form_pairwise_bundles");

    //> Copy the valid-pair count from GPU to host
    cudacheck(cudaMemcpy(&result.valid_pairs, dev_valid_pair_count, sizeof(unsigned), cudaMemcpyDeviceToHost));
    profile_lap(profiler, TimerCategory::DataTransfer, "bundle valid-pair count D->H");

    //> Free the GPU memory
    cudaFree(dev_valid_pair_count);
    profile_lap(profiler, TimerCategory::Other, "bundle temp counter free");

    return true;
}
