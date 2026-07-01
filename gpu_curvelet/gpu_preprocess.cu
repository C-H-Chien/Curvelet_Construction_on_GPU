#include "gpu_preprocess.hpp"
#include "gpu_common.hpp"
#include "device/gpu_preprocess_kernels.cuh"
#include "timer.hpp"

#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>

#include <algorithm>
#include <cstdio>
#include <vector>

namespace {

struct SpatialIndexHost {
    float *dev_edges = nullptr;
    long long *dev_cell_keys = nullptr;
    int *dev_edge_order = nullptr;
    long long *dev_unique_cells = nullptr;
    int *dev_cell_counts = nullptr;
    int *dev_cell_starts = nullptr;
    int *dev_sorted_edge_ids = nullptr;
    int num_cells = 0;
};

inline void profile_lap(CategoryProfiler *profiler, TimerCategory cat, const char *detail)
{
    if (profiler != nullptr) {
        profiler->lap(cat, detail);
    }
}

bool build_spatial_index(
    const GPUPreprocessConfig &cfg,
    const float *host_to_edges,
    SpatialIndexHost &index,
    CategoryProfiler *profiler)
{
    const int threads_per_block = 256;

    cudacheck(cudaMalloc(&index.dev_edges, static_cast<size_t>(cfg.num_edges) * cfg.sz_edge_data * sizeof(float)));
    cudacheck(cudaMalloc(&index.dev_cell_keys, static_cast<size_t>(cfg.num_edges) * sizeof(long long)));
    cudacheck(cudaMalloc(&index.dev_edge_order, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    cudacheck(cudaMalloc(&index.dev_unique_cells, static_cast<size_t>(cfg.num_edges) * sizeof(long long)));
    cudacheck(cudaMalloc(&index.dev_cell_counts, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    cudacheck(cudaMalloc(&index.dev_cell_starts, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    cudacheck(cudaMalloc(&index.dev_sorted_edge_ids, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    profile_lap(profiler, TimerCategory::MemoryAlloc, "spatial index buffers");

    cudacheck(cudaMemcpy(
        index.dev_edges, host_to_edges,
        static_cast<size_t>(cfg.num_edges) * cfg.sz_edge_data * sizeof(float),
        cudaMemcpyHostToDevice));
    profile_lap(profiler, TimerCategory::DataTransfer, "edges H->D");

    compute_cell_keys_kernel<<<div_up(cfg.num_edges, threads_per_block), threads_per_block>>>(
        cfg.num_edges, index.dev_edges, index.dev_cell_keys);
    cudacheck(cudaGetLastError());
    cudacheck(cudaDeviceSynchronize());
    profile_lap(profiler, TimerCategory::Kernel, "compute_cell_keys_kernel");

    thrust::device_ptr<long long> d_keys(index.dev_cell_keys);
    thrust::device_ptr<int> d_order(index.dev_edge_order);
    thrust::sequence(d_order, d_order + cfg.num_edges);
    thrust::sort_by_key(d_keys, d_keys + cfg.num_edges, d_order);

    thrust::device_ptr<long long> d_unique(index.dev_unique_cells);
    thrust::device_ptr<int> d_counts(index.dev_cell_counts);
    thrust::device_ptr<int> d_starts(index.dev_cell_starts);

    const auto ends = thrust::reduce_by_key( d_keys, d_keys + cfg.num_edges, thrust::make_constant_iterator(1), d_unique, d_counts );
    index.num_cells = static_cast<int>(ends.first - d_unique);

    thrust::exclusive_scan(d_counts, d_counts + index.num_cells, d_starts);
    profile_lap(profiler, TimerCategory::Thrust, "sort + reduce_by_key + scan");

    cudacheck(cudaMemcpy(
        index.dev_sorted_edge_ids, index.dev_edge_order,
        static_cast<size_t>(cfg.num_edges) * sizeof(int),
        cudaMemcpyDeviceToDevice));
    profile_lap(profiler, TimerCategory::DataTransfer, "sorted edge ids D->D");

    return true;
}

void free_spatial_index(SpatialIndexHost &index)
{
    if (index.dev_edges != nullptr) {
        cudaFree(index.dev_edges);
    }
    cudaFree(index.dev_cell_keys);
    cudaFree(index.dev_edge_order);
    cudaFree(index.dev_unique_cells);
    cudaFree(index.dev_cell_counts);
    cudaFree(index.dev_cell_starts);
    cudaFree(index.dev_sorted_edge_ids);
    index = SpatialIndexHost{};
}

SpatialIndexDevice spatial_device_view(const SpatialIndexHost &index)
{
    SpatialIndexDevice view;
    view.dev_edges = index.dev_edges;
    view.dev_unique_cells = index.dev_unique_cells;
    view.dev_cell_starts = index.dev_cell_starts;
    view.dev_cell_counts = index.dev_cell_counts;
    view.dev_sorted_edge_ids = index.dev_sorted_edge_ids;
    view.num_cells = index.num_cells;
    return view;
}

int finalize_neighbor_offsets(int num_edges, int *dev_neighbor_counts, int *dev_neighbor_offsets)
{
    thrust::device_ptr<int> d_counts(dev_neighbor_counts);
    thrust::device_ptr<int> d_offsets(dev_neighbor_offsets);
    thrust::exclusive_scan(d_counts, d_counts + num_edges, d_offsets);

    int last_offset = 0;
    int last_count = 0;
    cudacheck(cudaMemcpy( &last_offset, dev_neighbor_offsets + (num_edges - 1), sizeof(int), cudaMemcpyDeviceToHost ));
    cudacheck(cudaMemcpy( &last_count,  dev_neighbor_counts + (num_edges - 1),  sizeof(int), cudaMemcpyDeviceToHost ));
    int total_pairs = last_offset + last_count;
    cudacheck(cudaMemcpy( dev_neighbor_offsets + num_edges, &total_pairs, sizeof(int), cudaMemcpyHostToDevice ));

    return total_pairs;
}

void copy_csr_graph_to_host(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexHost &index,
    int total_pairs,
    int *dev_neighbor_offsets,
    int *dev_neighbor_ids,
    float *dev_neighbor_dist2,
    GPUNeighborGraph &graph)
{
    graph.host_edges = new float[static_cast<size_t>(cfg.num_edges) * cfg.sz_edge_data];
    graph.host_neighbor_offsets = new int[static_cast<size_t>(cfg.num_edges + 1)];
    if (total_pairs > 0) {
        graph.host_neighbor_ids = new int[static_cast<size_t>(total_pairs)];
        graph.host_neighbor_dist2 = new float[static_cast<size_t>(total_pairs)];
    }

    cudacheck(cudaMemcpy(graph.host_edges,            index.dev_edges,      static_cast<size_t>(cfg.num_edges) * cfg.sz_edge_data * sizeof(float), cudaMemcpyDeviceToHost));
    cudacheck(cudaMemcpy(graph.host_neighbor_offsets, dev_neighbor_offsets, static_cast<size_t>(cfg.num_edges + 1) * sizeof(int),                  cudaMemcpyDeviceToHost));
    if (total_pairs > 0) {
        cudacheck(cudaMemcpy( graph.host_neighbor_ids,   dev_neighbor_ids,   static_cast<size_t>(total_pairs) * sizeof(int),   cudaMemcpyDeviceToHost));
        cudacheck(cudaMemcpy( graph.host_neighbor_dist2, dev_neighbor_dist2, static_cast<size_t>(total_pairs) * sizeof(float), cudaMemcpyDeviceToHost));
    }
}

void print_csr_graph_stats(
    const char *label,
    const GPUPreprocessConfig &cfg,
    int total_pairs,
    unsigned max_num_of_neighbors)
{
    printf("%sCSR neighbor graph: %d anchors, %d total pairs, max number of neighbor edges %u\n", label, cfg.num_edges, total_pairs, max_num_of_neighbors);
    printf("%sCSR memory (device): edges %.2f MB + offsets %.2f KB + ids %.2f MB + dist2 %.2f MB = %.2f MB\n",
           label,
           static_cast<double>(cfg.num_edges) * cfg.sz_edge_data * sizeof(float) / (1024.0 * 1024.0),
           static_cast<double>(cfg.num_edges + 1) * sizeof(int) / 1024.0,
           static_cast<double>(total_pairs) * sizeof(int) / (1024.0 * 1024.0),
           static_cast<double>(total_pairs) * sizeof(float) / (1024.0 * 1024.0),
           (static_cast<double>(cfg.num_edges) * cfg.sz_edge_data * sizeof(float)
            + static_cast<double>(cfg.num_edges + 1) * sizeof(int)
            + static_cast<double>(total_pairs) * sizeof(int)
            + static_cast<double>(total_pairs) * sizeof(float)) / (1024.0 * 1024.0));
}

//> =============================== Kernel Launchers for the Two Pass CSR Strategy ===============================
//> Count the number of neighbors per anchor edge (used by the two pass CSR strategy)
void launch_neighbor_count_kernel(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexDevice &spatial,
    float rad_sqr,
    int *dev_neighbor_counts,
    unsigned *dev_max_num_of_neighbors,
    int threads_per_block)
{
    if (threads_per_block <= 0 || threads_per_block > 1024) {
        fprintf(stderr, "launch_neighbor_count_kernel: invalid threads_per_block=%d (use 1..1024)\n", threads_per_block);
        return;
    }

    const int num_blocks = div_up(cfg.num_edges, threads_per_block);

    count_neighbors_kernel<<<num_blocks, threads_per_block>>>( cfg.num_edges, spatial, rad_sqr, cfg.neighbor_radius, dev_neighbor_counts, dev_max_num_of_neighbors );
    cudacheck(cudaGetLastError());
}

//> Fill the CSR neighbor graph (used by the two pass CSR strategy)
void launch_fill_csr_neighbors_kernel(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexDevice &spatial,
    float rad_sqr,
    const int *dev_neighbor_counts,
    const int *dev_neighbor_offsets,
    int *dev_neighbor_ids,
    float *dev_neighbor_dist2,
    unsigned max_num_of_neighbors,
    int threads_per_block)
{
    if (threads_per_block <= 0 || threads_per_block > 1024) {
        fprintf(stderr, "launch_fill_csr_neighbors_kernel: invalid threads_per_block=%d (min: 1, max: 1024)\n", threads_per_block);
        return;
    }

    const int num_blocks = div_up(cfg.num_edges, threads_per_block);
    const size_t shmem_bytes = static_cast<size_t>(threads_per_block)
        * static_cast<size_t>(max_num_of_neighbors) * sizeof(NeighborCandidate);

    printf("[two-pass fill] threads/block=%d, max_degree=%u, shmem=%.2f KB\n",
           threads_per_block, max_num_of_neighbors, static_cast<double>(shmem_bytes) / 1024.0);

    fill_csr_neighbors_fixed_kernel<<<num_blocks, threads_per_block, shmem_bytes>>>(
        cfg.num_edges, max_num_of_neighbors, spatial, rad_sqr, cfg.neighbor_radius,
        dev_neighbor_counts, dev_neighbor_offsets, dev_neighbor_ids, dev_neighbor_dist2);
    cudacheck(cudaGetLastError());
}
//> ========================================================================================================================

void launch_discover_fixed_row_warp_kernel(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexDevice &spatial,
    float rad_sqr,
    int max_candidates,
    int *dev_neighbor_counts,
    int *dev_neighbor_list,
    float *dev_neighbor_dist2_row,
    unsigned *dev_max_num_of_neighbors,
    unsigned *dev_truncated_anchors)
{
    const int warps_per_block = cfg.neighbor_warps_per_block;
    if (warps_per_block <= 0 || warps_per_block > 32) {
        fprintf(stderr,
                "launch_discover_fixed_row_warp_kernel: invalid neighbor_warps_per_block=%d (use 1..32)\n",
                warps_per_block);
        return;
    }

    const int threads_per_block = warps_per_block * 32;
    const int num_blocks = div_up(cfg.num_edges, warps_per_block);
    const size_t shmem_bytes = static_cast<size_t>(warps_per_block) * static_cast<size_t>(max_candidates)
        * sizeof(NeighborCandidate)
        + static_cast<size_t>(warps_per_block) * sizeof(int);

    printf("[fixed-row warp] warps/block=%d, threads/block=%d, shmem=%.2f KB\n",
           warps_per_block, threads_per_block, static_cast<double>(shmem_bytes) / 1024.0);

    discover_fixed_row_warp_kernel<<<num_blocks, threads_per_block, shmem_bytes>>>(
        cfg.num_edges, max_candidates, max_candidates, warps_per_block,
        spatial, rad_sqr, cfg.neighbor_radius,
        dev_neighbor_counts, dev_neighbor_list, dev_neighbor_dist2_row,
        dev_max_num_of_neighbors, dev_truncated_anchors);
    cudacheck(cudaGetLastError());
}

void launch_discover_and_stage_neighbors_kernel(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexDevice &spatial,
    float rad_sqr,
    int max_candidates,
    int *dev_staged_counts,
    int *dev_staged_ids,
    float *dev_staged_dist2,
    unsigned *dev_max_num_of_neighbors,
    unsigned *dev_truncated_anchors,
    int threads_per_block)
{
    if (threads_per_block <= 0 || threads_per_block > 1024) {
        fprintf(stderr,
                "launch_discover_and_stage_neighbors_kernel: invalid threads_per_block=%d (use 1..1024)\n",
                threads_per_block);
        return;
    }

    const int num_blocks = div_up(cfg.num_edges, threads_per_block);
    const size_t shmem_bytes = static_cast<size_t>(threads_per_block)
        * static_cast<size_t>(max_candidates) * sizeof(NeighborCandidate);

    printf("[single-pass stage] threads/block=%d, shmem=%.2f KB\n",
           threads_per_block, static_cast<double>(shmem_bytes) / 1024.0);

    discover_and_stage_neighbors_kernel<<<num_blocks, threads_per_block, shmem_bytes>>>(
        cfg.num_edges, max_candidates, spatial, rad_sqr, cfg.neighbor_radius,
        dev_staged_counts, dev_staged_ids, dev_staged_dist2,
        dev_max_num_of_neighbors, dev_truncated_anchors);
    cudacheck(cudaGetLastError());
}

void launch_compact_staged_neighbors_kernel(
    const GPUPreprocessConfig &cfg,
    int max_candidates,
    const int *dev_staged_counts,
    const int *dev_staged_ids,
    const float *dev_staged_dist2,
    const int *dev_neighbor_offsets,
    int *dev_neighbor_ids,
    float *dev_neighbor_dist2,
    int threads_per_block)
{
    if (threads_per_block <= 0 || threads_per_block > 1024) {
        fprintf(stderr,
                "launch_compact_staged_neighbors_kernel: invalid threads_per_block=%d (use 1..1024)\n",
                threads_per_block);
        return;
    }

    const int num_blocks = div_up(cfg.num_edges, threads_per_block);

    compact_staged_neighbors_kernel<<<num_blocks, threads_per_block>>>(
        cfg.num_edges, max_candidates,
        dev_staged_counts, dev_staged_ids, dev_staged_dist2,
        dev_neighbor_offsets, dev_neighbor_ids, dev_neighbor_dist2);
    cudacheck(cudaGetLastError());
}

bool build_CSR_graph_twopass(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexHost &index,
    float rad_sqr,
    GPUNeighborGraph &graph,
    CategoryProfiler *profiler)
{
    const SpatialIndexDevice spatial = spatial_device_view(index);

    int *dev_neighbor_counts = nullptr;
    unsigned *dev_max_num_of_neighbors = nullptr;
    int *dev_neighbor_offsets = nullptr;
    cudacheck(cudaMalloc(&dev_neighbor_counts, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    cudacheck(cudaMalloc(&dev_max_num_of_neighbors, sizeof(unsigned)));
    cudacheck(cudaMalloc(&dev_neighbor_offsets, static_cast<size_t>(cfg.num_edges + 1) * sizeof(int)));
    cudacheck(cudaMemset(dev_max_num_of_neighbors, 0, sizeof(unsigned)));
    profile_lap(profiler, TimerCategory::MemoryAlloc, "two-pass count/offset buffers");

    launch_neighbor_count_kernel(
        cfg, spatial, rad_sqr,
        dev_neighbor_counts, dev_max_num_of_neighbors, cfg.neighbor_count_threads);
    cudacheck(cudaDeviceSynchronize());
    profile_lap(profiler, TimerCategory::Kernel, "two-pass count_neighbors_kernel");

    unsigned max_num_of_neighbors = 0;
    cudacheck(cudaMemcpy(&max_num_of_neighbors, dev_max_num_of_neighbors, sizeof(unsigned), cudaMemcpyDeviceToHost));
    profile_lap(profiler, TimerCategory::DataTransfer, "two-pass max degree D->H");

    const int total_pairs = finalize_neighbor_offsets(cfg.num_edges, dev_neighbor_counts, dev_neighbor_offsets);
    profile_lap(profiler, TimerCategory::Thrust, "two-pass neighbor offsets scan");

    int *dev_neighbor_ids = nullptr;
    float *dev_neighbor_dist2 = nullptr;
    if (total_pairs > 0) {
        cudacheck(cudaMalloc(&dev_neighbor_ids, static_cast<size_t>(total_pairs) * sizeof(int)));
        cudacheck(cudaMalloc(&dev_neighbor_dist2, static_cast<size_t>(total_pairs) * sizeof(float)));
    }
    profile_lap(profiler, TimerCategory::MemoryAlloc, "two-pass compact CSR buffers");

    launch_fill_csr_neighbors_kernel(
        cfg, spatial, rad_sqr, dev_neighbor_counts, dev_neighbor_offsets,
        dev_neighbor_ids, dev_neighbor_dist2, max_num_of_neighbors, cfg.neighbor_fill_threads);
    cudacheck(cudaDeviceSynchronize());
    profile_lap(profiler, TimerCategory::Kernel, "two-pass fill_csr_neighbors_kernel");

    graph.num_edges = cfg.num_edges;
    graph.total_neighbor_pairs = total_pairs;
    graph.max_neighbor_degree = max_num_of_neighbors;
    graph.layout = GPUNeighborLayout::CSR;
    graph.dev_edges = index.dev_edges;
    graph.dev_neighbor_offsets = dev_neighbor_offsets;
    graph.dev_neighbor_ids = dev_neighbor_ids;
    graph.dev_neighbor_dist2 = dev_neighbor_dist2;

    if (cfg.copy_to_host) {
        copy_csr_graph_to_host(
            cfg, index, total_pairs,
            dev_neighbor_offsets, dev_neighbor_ids, dev_neighbor_dist2, graph);
        profile_lap(profiler, TimerCategory::DataTransfer, "two-pass CSR D->H");
    }

    cudaFree(dev_neighbor_counts);
    cudaFree(dev_max_num_of_neighbors);
    profile_lap(profiler, TimerCategory::Other, "two-pass temp buffer free");

    print_csr_graph_stats("[two-pass] ", cfg, total_pairs, max_num_of_neighbors);
    return true;
}

bool build_CSR_graph_singlepass(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexHost &index,
    float rad_sqr,
    GPUNeighborGraph &graph,
    CategoryProfiler *profiler)
{
    const int max_candidates = static_cast<int>(cfg.max_candidates);
    const SpatialIndexDevice spatial = spatial_device_view(index);

    int *dev_staged_counts = nullptr;
    int *dev_staged_ids = nullptr;
    float *dev_staged_dist2 = nullptr;
    unsigned *dev_max_num_of_neighbors = nullptr;
    unsigned *dev_truncated_anchors = nullptr;
    int *dev_neighbor_offsets = nullptr;

    cudacheck(cudaMalloc(&dev_staged_counts, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    cudacheck(cudaMalloc(&dev_staged_ids,    static_cast<size_t>(cfg.num_edges) * static_cast<size_t>(max_candidates) * sizeof(int)));
    cudacheck(cudaMalloc(&dev_staged_dist2,  static_cast<size_t>(cfg.num_edges) * static_cast<size_t>(max_candidates) * sizeof(float)));
    cudacheck(cudaMalloc(&dev_max_num_of_neighbors, sizeof(unsigned)));
    cudacheck(cudaMalloc(&dev_truncated_anchors, sizeof(unsigned)));
    cudacheck(cudaMalloc(&dev_neighbor_offsets, static_cast<size_t>(cfg.num_edges + 1) * sizeof(int)));
    cudacheck(cudaMemset(dev_max_num_of_neighbors, 0, sizeof(unsigned)));
    cudacheck(cudaMemset(dev_truncated_anchors, 0, sizeof(unsigned)));
    profile_lap(profiler, TimerCategory::MemoryAlloc, "single-pass staging buffers");

    launch_discover_and_stage_neighbors_kernel(
        cfg, spatial, rad_sqr, max_candidates,
        dev_staged_counts, dev_staged_ids, dev_staged_dist2,
        dev_max_num_of_neighbors, dev_truncated_anchors, cfg.neighbor_stage_threads);
    cudacheck(cudaDeviceSynchronize());
    profile_lap(profiler, TimerCategory::Kernel, "single-pass discover_and_stage_neighbors_kernel");

    const int total_pairs = finalize_neighbor_offsets(cfg.num_edges, dev_staged_counts, dev_neighbor_offsets);
    profile_lap(profiler, TimerCategory::Thrust, "single-pass neighbor offsets scan");

    int *dev_neighbor_ids = nullptr;
    float *dev_neighbor_dist2 = nullptr;
    if (total_pairs > 0) {
        cudacheck(cudaMalloc(&dev_neighbor_ids, static_cast<size_t>(total_pairs) * sizeof(int)));
        cudacheck(cudaMalloc(&dev_neighbor_dist2, static_cast<size_t>(total_pairs) * sizeof(float)));
    }
    profile_lap(profiler, TimerCategory::MemoryAlloc, "single-pass compact CSR buffers");

    launch_compact_staged_neighbors_kernel(
        cfg, max_candidates,
        dev_staged_counts, dev_staged_ids, dev_staged_dist2,
        dev_neighbor_offsets, dev_neighbor_ids, dev_neighbor_dist2,
        cfg.neighbor_stage_threads);
    cudacheck(cudaDeviceSynchronize());
    profile_lap(profiler, TimerCategory::Kernel, "single-pass compact_staged_neighbors_kernel");

    unsigned max_num_of_neighbors = 0;
    unsigned truncated_anchors = 0;
    cudacheck(cudaMemcpy(&max_num_of_neighbors, dev_max_num_of_neighbors, sizeof(unsigned), cudaMemcpyDeviceToHost));
    cudacheck(cudaMemcpy(&truncated_anchors, dev_truncated_anchors, sizeof(unsigned), cudaMemcpyDeviceToHost));
    if (truncated_anchors > 0) {
        fprintf(stderr, "Warning: %u anchors hit max_candidates=%d; neighbor lists may be truncated\n",
                truncated_anchors, max_candidates);
    }
    profile_lap(profiler, TimerCategory::DataTransfer, "single-pass stats D->H");

    graph.num_edges = cfg.num_edges;
    graph.total_neighbor_pairs = total_pairs;
    graph.max_neighbor_degree = max_num_of_neighbors;
    graph.layout = GPUNeighborLayout::CSR;
    graph.dev_edges = index.dev_edges;
    graph.dev_neighbor_offsets = dev_neighbor_offsets;
    graph.dev_neighbor_ids = dev_neighbor_ids;
    graph.dev_neighbor_dist2 = dev_neighbor_dist2;

    if (cfg.copy_to_host) {
        copy_csr_graph_to_host(
            cfg, index, total_pairs,
            dev_neighbor_offsets, dev_neighbor_ids, dev_neighbor_dist2, graph);
        profile_lap(profiler, TimerCategory::DataTransfer, "single-pass CSR D->H");
    }

    cudaFree(dev_staged_counts);
    cudaFree(dev_staged_ids);
    cudaFree(dev_staged_dist2);
    cudaFree(dev_max_num_of_neighbors);
    cudaFree(dev_truncated_anchors);
    profile_lap(profiler, TimerCategory::Other, "single-pass staging buffer free");

    const double staged_mb = static_cast<double>(cfg.num_edges) * max_candidates * (sizeof(int) + sizeof(float)) / (1024.0 * 1024.0);
    printf("[single-pass] staged buffer: %.2f MB (num_edges x max_candidates)\n", staged_mb);
    print_csr_graph_stats("[single-pass] ", cfg, total_pairs, max_num_of_neighbors);
    return true;
}

void copy_fixed_row_graph_to_host(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexHost &index,
    int *dev_neighbor_counts,
    int *dev_neighbor_list,
    float *dev_neighbor_dist2_row,
    GPUNeighborGraph &graph)
{
    const size_t edge_bytes = static_cast<size_t>(cfg.num_edges) * cfg.sz_edge_data * sizeof(float);
    const int slots = graph.neighbor_slots_per_anchor;
    const size_t row_ints = static_cast<size_t>(cfg.num_edges) * static_cast<size_t>(slots);

    graph.host_edges = new float[static_cast<size_t>(cfg.num_edges) * cfg.sz_edge_data];
    graph.host_neighbor_counts = new int[static_cast<size_t>(cfg.num_edges)];
    graph.host_neighbor_list = new int[row_ints];
    graph.host_neighbor_dist2_row = new float[row_ints];

    cudacheck(cudaMemcpy(graph.host_edges, index.dev_edges, edge_bytes, cudaMemcpyDeviceToHost));
    cudacheck(cudaMemcpy(
        graph.host_neighbor_counts, dev_neighbor_counts,
        static_cast<size_t>(cfg.num_edges) * sizeof(int), cudaMemcpyDeviceToHost));
    cudacheck(cudaMemcpy(graph.host_neighbor_list, dev_neighbor_list, row_ints * sizeof(int), cudaMemcpyDeviceToHost));
    cudacheck(cudaMemcpy(
        graph.host_neighbor_dist2_row, dev_neighbor_dist2_row,
        row_ints * sizeof(float), cudaMemcpyDeviceToHost));
}

void print_fixed_row_graph_stats(
    const char *label,
    const GPUPreprocessConfig &cfg,
    int total_pairs,
    unsigned max_num_of_neighbors,
    int slots_per_anchor)
{
    const size_t row_elems = static_cast<size_t>(cfg.num_edges) * static_cast<size_t>(slots_per_anchor);
    printf("%sFixed-row neighbor graph: %d anchors, %d slots/anchor, %d total pairs, max degree %u\n",
           label, cfg.num_edges, slots_per_anchor, total_pairs, max_num_of_neighbors);
    printf("%sFixed-row memory (device): edges %.2f MB + list %.2f MB + dist2 %.2f MB + counts %.2f KB\n",
           label,
           static_cast<double>(cfg.num_edges) * cfg.sz_edge_data * sizeof(float) / (1024.0 * 1024.0),
           static_cast<double>(row_elems) * sizeof(int) / (1024.0 * 1024.0),
           static_cast<double>(row_elems) * sizeof(float) / (1024.0 * 1024.0),
           static_cast<double>(cfg.num_edges) * sizeof(int) / 1024.0);
}

bool build_neighbor_graph_fixed_row(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexHost &index,
    float rad_sqr,
    GPUNeighborGraph &graph,
    CategoryProfiler *profiler)
{
    const int max_candidates = static_cast<int>(cfg.max_candidates);
    const SpatialIndexDevice spatial = spatial_device_view(index);
    const size_t row_ints = static_cast<size_t>(cfg.num_edges) * static_cast<size_t>(max_candidates);

    int *dev_neighbor_counts = nullptr;
    int *dev_neighbor_list = nullptr;
    float *dev_neighbor_dist2_row = nullptr;
    unsigned *dev_max_num_of_neighbors = nullptr;
    unsigned *dev_truncated_anchors = nullptr;

    cudacheck(cudaMalloc(&dev_neighbor_counts, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    cudacheck(cudaMalloc(&dev_neighbor_list, row_ints * sizeof(int)));
    cudacheck(cudaMalloc(&dev_neighbor_dist2_row, row_ints * sizeof(float)));
    cudacheck(cudaMalloc(&dev_max_num_of_neighbors, sizeof(unsigned)));
    cudacheck(cudaMalloc(&dev_truncated_anchors, sizeof(unsigned)));
    cudacheck(cudaMemset(dev_neighbor_list, 0xFF, row_ints * sizeof(int)));
    cudacheck(cudaMemset(dev_neighbor_dist2_row, 0, row_ints * sizeof(float)));
    cudacheck(cudaMemset(dev_max_num_of_neighbors, 0, sizeof(unsigned)));
    cudacheck(cudaMemset(dev_truncated_anchors, 0, sizeof(unsigned)));
    profile_lap(profiler, TimerCategory::MemoryAlloc, "fixed-row stage neighbor buffers");

    launch_discover_and_stage_neighbors_kernel(
        cfg, spatial, rad_sqr, max_candidates,
        dev_neighbor_counts, dev_neighbor_list, dev_neighbor_dist2_row,
        dev_max_num_of_neighbors, dev_truncated_anchors, cfg.neighbor_stage_threads);
    cudacheck(cudaDeviceSynchronize());
    profile_lap(profiler, TimerCategory::Kernel, "fixed-row stage discover_and_stage_neighbors_kernel");

    unsigned max_num_of_neighbors = 0;
    unsigned truncated_anchors = 0;
    cudacheck(cudaMemcpy(&max_num_of_neighbors, dev_max_num_of_neighbors, sizeof(unsigned), cudaMemcpyDeviceToHost));
    cudacheck(cudaMemcpy(&truncated_anchors, dev_truncated_anchors, sizeof(unsigned), cudaMemcpyDeviceToHost));
    if (truncated_anchors > 0) {
        fprintf(stderr,
                "Warning: %u anchors hit max_candidates=%d; neighbor lists may be truncated\n",
                truncated_anchors, max_candidates);
    }
    profile_lap(profiler, TimerCategory::DataTransfer, "fixed-row stage stats D->H");

    int total_pairs = 0;
    if (cfg.copy_to_host) {
        std::vector<int> host_counts(static_cast<size_t>(cfg.num_edges));
        cudacheck(cudaMemcpy(
            host_counts.data(), dev_neighbor_counts,
            static_cast<size_t>(cfg.num_edges) * sizeof(int), cudaMemcpyDeviceToHost));
        for (int c : host_counts) {
            total_pairs += c;
        }
        profile_lap(profiler, TimerCategory::DataTransfer, "fixed-row stage counts D->H (host sum)");
    }
    else {
        thrust::device_ptr<int> d_counts(dev_neighbor_counts);
        total_pairs = thrust::reduce(d_counts, d_counts + cfg.num_edges, 0);
        profile_lap(profiler, TimerCategory::Thrust, "fixed-row stage reduce neighbor counts");
    }

    graph.layout = GPUNeighborLayout::FixedRow;
    graph.num_edges = cfg.num_edges;
    graph.total_neighbor_pairs = total_pairs;
    graph.max_neighbor_degree = max_num_of_neighbors;
    graph.neighbor_slots_per_anchor = max_candidates;
    graph.dev_edges = index.dev_edges;
    graph.dev_neighbor_list = dev_neighbor_list;
    graph.dev_neighbor_dist2_row = dev_neighbor_dist2_row;
    graph.dev_neighbor_counts = dev_neighbor_counts;

    if (cfg.copy_to_host) {
        copy_fixed_row_graph_to_host(
            cfg, index, dev_neighbor_counts, dev_neighbor_list, dev_neighbor_dist2_row, graph);
        profile_lap(profiler, TimerCategory::DataTransfer, "fixed-row stage graph D->H");
    }

    cudaFree(dev_max_num_of_neighbors);
    cudaFree(dev_truncated_anchors);
    profile_lap(profiler, TimerCategory::Other, "fixed-row stage temp buffer free");

    print_fixed_row_graph_stats("[fixed-row] ", cfg, total_pairs, max_num_of_neighbors, max_candidates);
    return true;
}

bool build_neighbor_graph_fixed_row_warp(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexHost &index,
    float rad_sqr,
    GPUNeighborGraph &graph,
    CategoryProfiler *profiler)
{
    const int max_candidates = static_cast<int>(cfg.max_candidates);
    const SpatialIndexDevice spatial = spatial_device_view(index);
    const size_t row_ints = static_cast<size_t>(cfg.num_edges) * static_cast<size_t>(max_candidates);

    int *dev_neighbor_counts = nullptr;
    int *dev_neighbor_list = nullptr;
    float *dev_neighbor_dist2_row = nullptr;
    unsigned *dev_max_num_of_neighbors = nullptr;
    unsigned *dev_truncated_anchors = nullptr;

    cudacheck(cudaMalloc(&dev_neighbor_counts, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    cudacheck(cudaMalloc(&dev_neighbor_list, row_ints * sizeof(int)));
    cudacheck(cudaMalloc(&dev_neighbor_dist2_row, row_ints * sizeof(float)));
    cudacheck(cudaMalloc(&dev_max_num_of_neighbors, sizeof(unsigned)));
    cudacheck(cudaMalloc(&dev_truncated_anchors, sizeof(unsigned)));
    cudacheck(cudaMemset(dev_neighbor_list, 0xFF, row_ints * sizeof(int)));
    cudacheck(cudaMemset(dev_neighbor_dist2_row, 0, row_ints * sizeof(float)));
    cudacheck(cudaMemset(dev_max_num_of_neighbors, 0, sizeof(unsigned)));
    cudacheck(cudaMemset(dev_truncated_anchors, 0, sizeof(unsigned)));
    profile_lap(profiler, TimerCategory::MemoryAlloc, "fixed-row warp neighbor buffers");

    launch_discover_fixed_row_warp_kernel(
        cfg, spatial, rad_sqr, max_candidates,
        dev_neighbor_counts, dev_neighbor_list, dev_neighbor_dist2_row,
        dev_max_num_of_neighbors, dev_truncated_anchors);
    cudacheck(cudaDeviceSynchronize());
    profile_lap(profiler, TimerCategory::Kernel, "fixed-row warp discover_fixed_row_warp_kernel");

    unsigned max_num_of_neighbors = 0;
    unsigned truncated_anchors = 0;
    cudacheck(cudaMemcpy(&max_num_of_neighbors, dev_max_num_of_neighbors, sizeof(unsigned), cudaMemcpyDeviceToHost));
    cudacheck(cudaMemcpy(&truncated_anchors, dev_truncated_anchors, sizeof(unsigned), cudaMemcpyDeviceToHost));
    if (truncated_anchors > 0) {
        fprintf(stderr,
                "Warning: %u anchors hit max_candidates=%d; neighbor lists may be truncated\n",
                truncated_anchors, max_candidates);
    }
    profile_lap(profiler, TimerCategory::DataTransfer, "fixed-row warp stats D->H");

    int total_pairs = 0;
    if (cfg.copy_to_host) {
        std::vector<int> host_counts(static_cast<size_t>(cfg.num_edges));
        cudacheck(cudaMemcpy( host_counts.data(), dev_neighbor_counts, static_cast<size_t>(cfg.num_edges) * sizeof(int), cudaMemcpyDeviceToHost ));
        for (int c : host_counts) {
            total_pairs += c;
        }
        profile_lap(profiler, TimerCategory::DataTransfer, "fixed-row warp counts D->H (host sum)");
    }
    else {
        thrust::device_ptr<int> d_counts(dev_neighbor_counts);
        total_pairs = thrust::reduce(d_counts, d_counts + cfg.num_edges, 0);
        profile_lap(profiler, TimerCategory::Thrust, "fixed-row warp reduce neighbor counts");
    }

    graph.layout = GPUNeighborLayout::FixedRow;
    graph.num_edges = cfg.num_edges;
    graph.total_neighbor_pairs = total_pairs;
    graph.max_neighbor_degree = max_num_of_neighbors;
    graph.neighbor_slots_per_anchor = max_candidates;
    graph.dev_edges = index.dev_edges;
    graph.dev_neighbor_list = dev_neighbor_list;
    graph.dev_neighbor_dist2_row = dev_neighbor_dist2_row;
    graph.dev_neighbor_counts = dev_neighbor_counts;

    if (cfg.copy_to_host) {
        copy_fixed_row_graph_to_host(cfg, index, dev_neighbor_counts, dev_neighbor_list, dev_neighbor_dist2_row, graph);
        profile_lap(profiler, TimerCategory::DataTransfer, "fixed-row warp graph D->H");
    }

    cudaFree(dev_max_num_of_neighbors);
    cudaFree(dev_truncated_anchors);
    profile_lap(profiler, TimerCategory::Other, "fixed-row warp temp buffer free");

    print_fixed_row_graph_stats("[fixed-row warp] ", cfg, total_pairs, max_num_of_neighbors, max_candidates);
    return true;
}

bool build_CSR_graph(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexHost &index,
    float rad_sqr,
    GPUNeighborGraph &graph,
    CategoryProfiler *profiler)
{
    if (cfg.csr_strategy == GPUNeighborCSRStrategy::TwoPass) {
        std::cout << "Building CSR graph using two-pass strategy" << std::endl;
        return build_CSR_graph_twopass(cfg, index, rad_sqr, graph, profiler);
    }

    return build_CSR_graph_singlepass(cfg, index, rad_sqr, graph, profiler);
}

bool build_neighbor_graph(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexHost &index,
    float rad_sqr,
    GPUNeighborGraph &graph,
    CategoryProfiler *profiler)
{
    if (cfg.neighbor_layout == GPUNeighborLayout::FixedRow) {
        if (cfg.fixed_row_build == GPUFixedRowBuildStrategy::Warp) {
            std::cout << "Building neighbor graph using fixed-row warp layout (one warp per anchor)" << std::endl;
            return build_neighbor_graph_fixed_row_warp(cfg, index, rad_sqr, graph, profiler);
        }
        std::cout << "Building neighbor graph using fixed-row layout (edgeLookList-style)" << std::endl;
        return build_neighbor_graph_fixed_row(cfg, index, rad_sqr, graph, profiler);
    }

    return build_CSR_graph(cfg, index, rad_sqr, graph, profiler);
}

} // namespace

//> MAIN: Build the CSR graph for the neighbor-search stage of curvelet construction.
bool gpu_preprocess_build(
    const GPUPreprocessConfig &cfg,
    const float *host_to_edges,
    GPUPreprocessResult &result,
    CategoryProfiler *profiler)
{
    //> Free the previous result just in case
    gpu_preprocess_free(result);

    if (cfg.num_edges <= 0 || cfg.sz_edge_data < kEdgeFields || host_to_edges == nullptr) {
        fprintf(stderr,
                "gpu_preprocess_build: invalid inputs: num_edges=%d, sz_edge_data=%d, host_to_edges=%p\n",
                cfg.num_edges, cfg.sz_edge_data, static_cast<const void *>(host_to_edges));
        return false;
    }

    cudacheck(cudaSetDevice(cfg.device_id));

    //> Radius centered at anchor edge
    const float rad_sqr = cfg.rad * cfg.rad;

    //> PHASE I: BUILD SPATIAL INDEX
    SpatialIndexHost spatial;
    if (!build_spatial_index(cfg, host_to_edges, spatial, profiler)) {
        return false;
    }

    //> PHASE II: BUILD NEIGHBOR GRAPH
    if (!build_neighbor_graph(cfg, spatial, rad_sqr, result.csr, profiler)) {
        free_spatial_index(spatial);
        return false;
    }

    //> dev_edges ownership transferred to result.csr; clear pointer so free_spatial_index skips it
    spatial.dev_edges = nullptr;
    free_spatial_index(spatial);
    profile_lap(profiler, TimerCategory::Other, "free_spatial_index");

    return true;
}

void gpu_preprocess_free(GPUNeighborGraph &graph)
{
    if (graph.dev_edges != nullptr) {
        cudaFree(graph.dev_edges);
        graph.dev_edges = nullptr;
    }

    cudaFree(graph.dev_neighbor_offsets);
    cudaFree(graph.dev_neighbor_ids);
    cudaFree(graph.dev_neighbor_dist2);
    cudaFree(graph.dev_neighbor_list);
    cudaFree(graph.dev_neighbor_dist2_row);
    cudaFree(graph.dev_neighbor_counts);

    delete[] graph.host_edges;
    delete[] graph.host_neighbor_offsets;
    delete[] graph.host_neighbor_ids;
    delete[] graph.host_neighbor_dist2;
    delete[] graph.host_neighbor_list;
    delete[] graph.host_neighbor_dist2_row;
    delete[] graph.host_neighbor_counts;
}

void gpu_preprocess_free(GPUPreprocessResult &result)
{
    gpu_preprocess_free(result.csr);
    result.csr = GPUNeighborGraph{};
}
