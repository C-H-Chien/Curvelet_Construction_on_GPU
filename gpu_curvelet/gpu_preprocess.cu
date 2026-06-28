#include "gpu_preprocess.hpp"
#include "gpu_common.hpp"
#include "device/gpu_preprocess_kernels.cuh"

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

bool build_spatial_index(
    const GPUPreprocessConfig &cfg,
    const float *host_to_edges,
    SpatialIndexHost &index)
{
    const int threads_per_block = 256;

    //> Allocate GPU memory
    cudacheck(cudaMalloc(&index.dev_edges, static_cast<size_t>(cfg.num_edges) * cfg.sz_edge_data * sizeof(float)));
    cudacheck(cudaMalloc(&index.dev_cell_keys, static_cast<size_t>(cfg.num_edges) * sizeof(long long)));
    cudacheck(cudaMalloc(&index.dev_edge_order, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    cudacheck(cudaMalloc(&index.dev_unique_cells, static_cast<size_t>(cfg.num_edges) * sizeof(long long)));
    cudacheck(cudaMalloc(&index.dev_cell_counts, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    cudacheck(cudaMalloc(&index.dev_cell_starts, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    cudacheck(cudaMalloc(&index.dev_sorted_edge_ids, static_cast<size_t>(cfg.num_edges) * sizeof(int)));

    cudacheck(cudaMemcpy(
        index.dev_edges, host_to_edges,
        static_cast<size_t>(cfg.num_edges) * cfg.sz_edge_data * sizeof(float),
        cudaMemcpyHostToDevice));

    //> Kernel 1: compute cell keys (subpixel edge locations to integer-based pixel cells)
    compute_cell_keys_kernel<<<div_up(cfg.num_edges, threads_per_block), threads_per_block>>>(
        cfg.num_edges, index.dev_edges, index.dev_cell_keys);
    cudacheck(cudaGetLastError());

    //> Kernel 2: sort the edges by cell keys
    //> Since the dev_cell_keys are packed 64-bit long long integers, we sort them so edges belonging the same cell become contiguous
    //> The cells appear in increasing order of the packed 64-bit integer key
    thrust::device_ptr<long long> d_keys(index.dev_cell_keys);
    thrust::device_ptr<int> d_order(index.dev_edge_order);
    thrust::sequence(d_order, d_order + cfg.num_edges);
    thrust::sort_by_key(d_keys, d_keys + cfg.num_edges, d_order);

    //> Then, we find the unique cells and count the number of edges in each unique cell, which gives the
    //  "offsets/starts" to the start of each cell
    //> Example:
    //  sorted keys:  [A, A, A, B, B, C, ...]
    //  unique:       [A,    B,    C]
    //  counts:       [3,    2,    1]
    //  starts:       [0,    3,    5]
    thrust::device_ptr<long long> d_unique(index.dev_unique_cells);
    thrust::device_ptr<int> d_counts(index.dev_cell_counts);
    thrust::device_ptr<int> d_starts(index.dev_cell_starts);

    //> Count the number of edges in each unique cell
    const auto ends = thrust::reduce_by_key(
        d_keys, d_keys + cfg.num_edges, thrust::make_constant_iterator(1), d_unique, d_counts);
    index.num_cells = static_cast<int>(ends.first - d_unique);

    //> Exclusive scan to get the starting index of each cell
    thrust::exclusive_scan(d_counts, d_counts + index.num_cells, d_starts);

    //> Copy the edge order (on device) to the sorted edge ids (which is also on device)
    //> TODO: Maybe we can stay with the same device memory to avoid memory copy?
    cudacheck(cudaMemcpy(
        index.dev_sorted_edge_ids, index.dev_edge_order,
        static_cast<size_t>(cfg.num_edges) * sizeof(int),
        cudaMemcpyDeviceToDevice));

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
    if (num_edges <= 0) {
        return 0;
    }

    thrust::device_ptr<int> d_counts(dev_neighbor_counts);
    thrust::device_ptr<int> d_offsets(dev_neighbor_offsets);
    thrust::exclusive_scan(d_counts, d_counts + num_edges, d_offsets);

    int total_pairs = 0;
    if (num_edges > 0) {
        int last_offset = 0;
        int last_count = 0;
        cudacheck(cudaMemcpy(
            &last_offset,
            dev_neighbor_offsets + (num_edges - 1),
            sizeof(int),
            cudaMemcpyDeviceToHost));
        cudacheck(cudaMemcpy(
            &last_count,
            dev_neighbor_counts + (num_edges - 1),
            sizeof(int),
            cudaMemcpyDeviceToHost));
        total_pairs = last_offset + last_count;
        cudacheck(cudaMemcpy(
            dev_neighbor_offsets + num_edges,
            &total_pairs,
            sizeof(int),
            cudaMemcpyHostToDevice));
    }
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
    const size_t edge_bytes = static_cast<size_t>(cfg.num_edges) * cfg.sz_edge_data * sizeof(float);
    graph.host_edges = new float[static_cast<size_t>(cfg.num_edges) * cfg.sz_edge_data];
    graph.host_neighbor_offsets = new int[static_cast<size_t>(cfg.num_edges + 1)];
    if (total_pairs > 0) {
        graph.host_neighbor_ids = new int[static_cast<size_t>(total_pairs)];
        graph.host_neighbor_dist2 = new float[static_cast<size_t>(total_pairs)];
    }

    cudacheck(cudaMemcpy(graph.host_edges, index.dev_edges, edge_bytes, cudaMemcpyDeviceToHost));
    cudacheck(cudaMemcpy(
        graph.host_neighbor_offsets,
        dev_neighbor_offsets,
        static_cast<size_t>(cfg.num_edges + 1) * sizeof(int),
        cudaMemcpyDeviceToHost));
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

bool build_CSR_graph_twopass(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexHost &index,
    float rad_sqr,
    GPUNeighborGraph &graph)
{
    const int max_candidates = static_cast<int>(cfg.max_candidates);
    const SpatialIndexDevice spatial = spatial_device_view(index);
    const size_t shmem_bytes = static_cast<size_t>(max_candidates) * sizeof(NeighborCandidate);

    int *dev_neighbor_counts = nullptr;
    unsigned *dev_max_num_of_neighbors = nullptr;
    int *dev_neighbor_offsets = nullptr;
    cudacheck(cudaMalloc(&dev_neighbor_counts, static_cast<size_t>(cfg.num_edges) * sizeof(int)));
    cudacheck(cudaMalloc(&dev_max_num_of_neighbors, sizeof(unsigned)));
    cudacheck(cudaMalloc(&dev_neighbor_offsets, static_cast<size_t>(cfg.num_edges + 1) * sizeof(int)));
    //> Initialize the maximal number of neighbors to 0
    cudacheck(cudaMemset(dev_max_num_of_neighbors, 0, sizeof(unsigned)));

    //> Kernel 3: identify neighbor edge cells and count the number of neighbors for each anchor edge
    count_neighbors_kernel<<<cfg.num_edges, 1, shmem_bytes>>>(
        cfg.num_edges, max_candidates, spatial, rad_sqr, cfg.neighbor_radius,
        dev_neighbor_counts, dev_max_num_of_neighbors);
    cudacheck(cudaGetLastError());

    //> Get the offsets of the number of neighbors in the `dev_neighbor_counts` array for all anchor edges.
    //> Also get the total number of anchor <-> neighbor pairs across all anchors; this is the size of
    //  `dev_neighbor_ids[]` and `dev_neighbor_dist2[]`.
    //> Copy the total number of pairs to the last element of the neighbor offsets. This enables accessing
    //  the neighbors of an anchor edge i from
    //  `dev_neighbor_ids[dev_neighbor_offsets[i] ... dev_neighbor_offsets[i+1])`, and
    //  `dev_neighbor_dist2[dev_neighbor_offsets[i] ... dev_neighbor_offsets[i+1])`.
    const int total_pairs = finalize_neighbor_offsets(
        cfg.num_edges, dev_neighbor_counts, dev_neighbor_offsets);

    //> Allocate memory for the neighbor ids and distances
    int *dev_neighbor_ids = nullptr;
    float *dev_neighbor_dist2 = nullptr;
    if (total_pairs > 0) {
        cudacheck(cudaMalloc(&dev_neighbor_ids, static_cast<size_t>(total_pairs) * sizeof(int)));
        cudacheck(cudaMalloc(&dev_neighbor_dist2, static_cast<size_t>(total_pairs) * sizeof(float)));
    }

    //> Kernel 4: fill the CSR neighbors for each anchor edge
    fill_csr_neighbors_kernel<<<cfg.num_edges, 1, shmem_bytes>>>(
        cfg.num_edges, max_candidates, spatial, rad_sqr, cfg.neighbor_radius,
        dev_neighbor_offsets, dev_neighbor_ids, dev_neighbor_dist2);
    cudacheck(cudaGetLastError());
    cudacheck(cudaDeviceSynchronize());

    //> Get the maximum amount of neighbors among all anchor edges
    unsigned max_num_of_neighbors = 0;
    cudacheck(cudaMemcpy(
        &max_num_of_neighbors, dev_max_num_of_neighbors, sizeof(unsigned), cudaMemcpyDeviceToHost));

    //> Fill in the graph structure
    graph.num_edges = cfg.num_edges;
    graph.total_neighbor_pairs = total_pairs;
    graph.max_neighbor_degree = max_num_of_neighbors;
    graph.dev_edges = index.dev_edges;
    graph.dev_neighbor_offsets = dev_neighbor_offsets;
    graph.dev_neighbor_ids = dev_neighbor_ids;
    graph.dev_neighbor_dist2 = dev_neighbor_dist2;

    //> This is meant for validation purposes
    if (cfg.copy_to_host) {
        copy_csr_graph_to_host(
            cfg, index, total_pairs,
            dev_neighbor_offsets, dev_neighbor_ids, dev_neighbor_dist2, graph);
    }

    cudaFree(dev_neighbor_counts);
    cudaFree(dev_max_num_of_neighbors);

    print_csr_graph_stats("[two-pass] ", cfg, total_pairs, max_num_of_neighbors);
    return true;
}

bool build_CSR_graph_singlepass(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexHost &index,
    float rad_sqr,
    GPUNeighborGraph &graph)
{
    const int max_candidates = static_cast<int>(cfg.max_candidates);
    const SpatialIndexDevice spatial = spatial_device_view(index);
    const size_t shmem_bytes = static_cast<size_t>(max_candidates) * sizeof(NeighborCandidate);

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

    discover_and_stage_neighbors_kernel<<<cfg.num_edges, 1, shmem_bytes>>>(
        cfg.num_edges, max_candidates, spatial, rad_sqr, cfg.neighbor_radius,
        dev_staged_counts, dev_staged_ids, dev_staged_dist2,
        dev_max_num_of_neighbors, dev_truncated_anchors);
    cudacheck(cudaGetLastError());

    const int total_pairs = finalize_neighbor_offsets( cfg.num_edges, dev_staged_counts, dev_neighbor_offsets );

    int *dev_neighbor_ids = nullptr;
    float *dev_neighbor_dist2 = nullptr;
    if (total_pairs > 0) {
        cudacheck(cudaMalloc(&dev_neighbor_ids, static_cast<size_t>(total_pairs) * sizeof(int)));
        cudacheck(cudaMalloc(&dev_neighbor_dist2, static_cast<size_t>(total_pairs) * sizeof(float)));
    }

    compact_staged_neighbors_kernel<<<cfg.num_edges, 1>>>(
        cfg.num_edges, max_candidates,
        dev_staged_counts, dev_staged_ids, dev_staged_dist2,
        dev_neighbor_offsets, dev_neighbor_ids, dev_neighbor_dist2);
    cudacheck(cudaGetLastError());
    cudacheck(cudaDeviceSynchronize());

    unsigned max_num_of_neighbors = 0;
    unsigned truncated_anchors = 0;
    cudacheck(cudaMemcpy( &max_num_of_neighbors, dev_max_num_of_neighbors, sizeof(unsigned), cudaMemcpyDeviceToHost) );
    cudacheck(cudaMemcpy( &truncated_anchors, dev_truncated_anchors, sizeof(unsigned), cudaMemcpyDeviceToHost) );
    if (truncated_anchors > 0) {
        fprintf(stderr, "Warning: %u anchors hit max_candidates=%d; neighbor lists may be truncated\n", truncated_anchors, max_candidates);
    }

    graph.num_edges = cfg.num_edges;
    graph.total_neighbor_pairs = total_pairs;
    graph.max_neighbor_degree = max_num_of_neighbors;
    graph.dev_edges = index.dev_edges;
    graph.dev_neighbor_offsets = dev_neighbor_offsets;
    graph.dev_neighbor_ids = dev_neighbor_ids;
    graph.dev_neighbor_dist2 = dev_neighbor_dist2;

    if (cfg.copy_to_host) {
        copy_csr_graph_to_host(
            cfg, index, total_pairs,
            dev_neighbor_offsets, dev_neighbor_ids, dev_neighbor_dist2, graph);
    }

    cudaFree(dev_staged_counts);
    cudaFree(dev_staged_ids);
    cudaFree(dev_staged_dist2);
    cudaFree(dev_max_num_of_neighbors);
    cudaFree(dev_truncated_anchors);

    const double staged_mb = static_cast<double>(cfg.num_edges) * max_candidates * (sizeof(int) + sizeof(float)) / (1024.0 * 1024.0);
    printf("[single-pass] staged buffer: %.2f MB (num_edges x max_candidates)\n", staged_mb);
    print_csr_graph_stats("[single-pass] ", cfg, total_pairs, max_num_of_neighbors);
    return true;
}

void free_csr_graph_device_only(GPUNeighborGraph &graph)
{
    cudaFree(graph.dev_neighbor_offsets);
    cudaFree(graph.dev_neighbor_ids);
    cudaFree(graph.dev_neighbor_dist2);
    delete[] graph.host_edges;
    delete[] graph.host_neighbor_offsets;
    delete[] graph.host_neighbor_ids;
    delete[] graph.host_neighbor_dist2;
    graph.dev_neighbor_offsets = nullptr;
    graph.dev_neighbor_ids = nullptr;
    graph.dev_neighbor_dist2 = nullptr;
    graph.host_edges = nullptr;
    graph.host_neighbor_offsets = nullptr;
    graph.host_neighbor_ids = nullptr;
    graph.host_neighbor_dist2 = nullptr;
}

bool build_CSR_graph(
    const GPUPreprocessConfig &cfg,
    const SpatialIndexHost &index,
    float rad_sqr,
    GPUNeighborGraph &graph)
{
    if (cfg.csr_strategy == GPUNeighborCSRStrategy::TwoPass) {
        std::cout << "Building CSR graph using two-pass strategy" << std::endl;
        return build_CSR_graph_twopass(cfg, index, rad_sqr, graph);
    }

    if (cfg.csr_strategy == GPUNeighborCSRStrategy::CompareCSR) {
        GPUPreprocessConfig single_cfg = cfg;
        GPUPreprocessConfig two_cfg = cfg;
        single_cfg.csr_strategy = GPUNeighborCSRStrategy::SinglePass;
        two_cfg.csr_strategy = GPUNeighborCSRStrategy::TwoPass;

        if (!build_CSR_graph_singlepass(single_cfg, index, rad_sqr, graph)) {
            return false;
        }

        GPUNeighborGraph twopass{};
        twopass.dev_edges = index.dev_edges;
        if (!build_CSR_graph_twopass(two_cfg, index, rad_sqr, twopass)) {
            return false;
        }

        if (cfg.copy_to_host) {
            const bool match = gpu_preprocess_compare_csr_graphs(graph, twopass, cfg.num_edges);
            printf("single-pass vs two-pass CSR comparison: %s\n", match ? "MATCH" : "MISMATCH");
        }

        free_csr_graph_device_only(twopass);
        return true;
    }

    return build_CSR_graph_singlepass(cfg, index, rad_sqr, graph);
}

} // namespace

bool gpu_preprocess_build(
    const GPUPreprocessConfig &cfg,
    const float *host_to_edges,
    GPUPreprocessResult &result)
{
    gpu_preprocess_free(result);

    if (cfg.num_edges <= 0 || cfg.sz_edge_data < kEdgeFields || host_to_edges == nullptr) {
        fprintf(stderr,
                "gpu_preprocess_build: invalid inputs: num_edges=%d, sz_edge_data=%d, host_to_edges=%p\n",
                cfg.num_edges, cfg.sz_edge_data, static_cast<const void *>(host_to_edges));
        return false;
    }

    cudacheck(cudaSetDevice(cfg.device_id));

    const float rad_sqr = cfg.rad * cfg.rad;

    //> PHASE I: BUILD SPATIAL INDEX
    SpatialIndexHost spatial;
    if (!build_spatial_index(cfg, host_to_edges, spatial)) {
        return false;
    }

    //> PHASE II: BUILD CSR GRAPH
    if (!build_CSR_graph(cfg, spatial, rad_sqr, result.csr)) {
        free_spatial_index(spatial);
        return false;
    }

    //> dev_edges ownership transferred to result.csr; clear pointer so free_spatial_index skips it
    spatial.dev_edges = nullptr;
    free_spatial_index(spatial);

    return true;
}

void gpu_preprocess_free(GPUPreprocessResult &result)
{
    if (result.csr.dev_edges != nullptr) {
        cudaFree(result.csr.dev_edges);
        result.csr.dev_edges = nullptr;
    }
    cudaFree(result.csr.dev_neighbor_offsets);
    cudaFree(result.csr.dev_neighbor_ids);
    cudaFree(result.csr.dev_neighbor_dist2);
    delete[] result.csr.host_edges;
    delete[] result.csr.host_neighbor_offsets;
    delete[] result.csr.host_neighbor_ids;
    delete[] result.csr.host_neighbor_dist2;
    result.csr = GPUNeighborGraph{};
}

bool gpu_preprocess_compare_csr_graphs(
    const GPUNeighborGraph &a,
    const GPUNeighborGraph &b,
    int num_edges_to_check)
{
    if (a.host_neighbor_offsets == nullptr || a.host_neighbor_ids == nullptr
        || b.host_neighbor_offsets == nullptr || b.host_neighbor_ids == nullptr) {
        fprintf(stderr, "gpu_preprocess_compare_csr_graphs: need host copies of both graphs\n");
        return false;
    }

    const int num_edges = std::min(a.num_edges, b.num_edges);
    const int check_n = (num_edges_to_check < 0) ? num_edges : std::min(num_edges_to_check, num_edges);

    int mismatches = 0;
    for (int anchor = 0; anchor < check_n; ++anchor) {
        const int a_begin = a.host_neighbor_offsets[anchor];
        const int a_end = a.host_neighbor_offsets[anchor + 1];
        const int b_begin = b.host_neighbor_offsets[anchor];
        const int b_end = b.host_neighbor_offsets[anchor + 1];
        const int a_count = a_end - a_begin;
        const int b_count = b_end - b_begin;

        if (a_count != b_count) {
            if (mismatches < 5) {
                fprintf(stderr,
                        "  anchor %d: count A=%d B=%d\n",
                        anchor, a_count, b_count);
            }
            ++mismatches;
            continue;
        }

        for (int k = 0; k < a_count; ++k) {
            const int a_id = a.host_neighbor_ids[a_begin + k];
            const int b_id = b.host_neighbor_ids[b_begin + k];
            if (a_id != b_id) {
                if (mismatches < 5) {
                    fprintf(stderr,
                            "  anchor %d neighbor %d: ids %d vs %d\n",
                            anchor, k, a_id, b_id);
                }
                ++mismatches;
            }
        }
    }

    if (mismatches == 0) {
        //> All neighbor ids and counts match
        printf("  compared %d anchors: CSR graphs match\n", check_n);
    } 
    else {
        //> Found mismatches across some anchors
        fprintf(stderr, "  found %d CSR mismatches across %d anchors\n", mismatches, check_n);
    }
    return mismatches == 0;
}
