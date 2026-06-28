#ifndef GPU_PREPROCESS_KERNELS_CUH
#define GPU_PREPROCESS_KERNELS_CUH

#include "gpu_dev_neighbors.cuh"

__global__ void compute_cell_keys_kernel(
    int num_edges,
    const float *dev_edges,
    long long *dev_cell_keys)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_edges) {
        return;
    }
    const float x = dev_edges[i * kEdgeFields + 0];
    const float y = dev_edges[i * kEdgeFields + 1];
    const int cx = static_cast<int>(lrintf(x));
    const int cy = static_cast<int>(lrintf(y));

    //> Stack cx and cy as a 64-bit integer
    dev_cell_keys[i] = pack_cell_key(cx, cy);
}

__global__ void count_neighbors_kernel(
    int num_edges,
    int max_candidates,
    SpatialIndexDevice spatial,
    float rad_sqr,
    unsigned neighbor_radius,
    int *dev_neighbor_counts,
    unsigned *dev_max_num_of_neighbors)
{
    extern __shared__ NeighborCandidate s_candidates[];

    const int te_idx = blockIdx.x;
    if (te_idx >= num_edges) {
        return;
    }

    const int count = discover_neighbors_by_binary_search(
        te_idx, spatial.dev_edges, spatial,
        rad_sqr, neighbor_radius,
        s_candidates, max_candidates);

    dev_neighbor_counts[te_idx] = count;
    if (count > 0) {
        atomicMax(dev_max_num_of_neighbors, static_cast<unsigned>(count));
    }
}

__global__ void fill_csr_neighbors_kernel(
    int num_edges,
    int max_candidates,
    SpatialIndexDevice spatial,
    float rad_sqr,
    unsigned neighbor_radius,
    const int *dev_neighbor_offsets,
    int *dev_neighbor_ids,
    float *dev_neighbor_dist2)
{
    extern __shared__ NeighborCandidate s_candidates[];

    const int te_idx = blockIdx.x;
    if (te_idx >= num_edges) {
        return;
    }

    int count = discover_neighbors_by_binary_search(
        te_idx, spatial.dev_edges, spatial, rad_sqr, neighbor_radius,
        s_candidates, max_candidates);
    sort_neighbors_by_distance(s_candidates, count);

    const int base = dev_neighbor_offsets[te_idx];
    for (int i = 0; i < count; ++i) {
        dev_neighbor_ids[base + i] = s_candidates[i].edge_id;
        dev_neighbor_dist2[base + i] = s_candidates[i].dist;
    }
}

//> Single-pass CSR: discover + sort once, stage ids/dist in a fixed per-anchor buffer.
__global__ void discover_and_stage_neighbors_kernel(
    int num_edges,
    int max_candidates,
    SpatialIndexDevice spatial,
    float rad_sqr,
    unsigned neighbor_radius,
    int *dev_staged_counts,
    int *dev_staged_ids,
    float *dev_staged_dist2,
    unsigned *dev_max_num_of_neighbors,
    unsigned *dev_truncated_anchors)
{
    extern __shared__ NeighborCandidate s_candidates[];

    const int te_idx = blockIdx.x;
    if (te_idx >= num_edges) {
        return;
    }

    const int raw_count = discover_neighbors_by_binary_search(
        te_idx, spatial.dev_edges, spatial, rad_sqr, neighbor_radius,
        s_candidates, max_candidates);
    sort_neighbors_by_distance(s_candidates, raw_count);

    int count = raw_count;
    if (raw_count >= max_candidates) {
        count = max_candidates;
        atomicAdd(dev_truncated_anchors, 1u);
    }

    dev_staged_counts[te_idx] = count;
    if (count > 0) {
        atomicMax(dev_max_num_of_neighbors, static_cast<unsigned>(count));
    }

    const int base = te_idx * max_candidates;
    for (int i = 0; i < count; ++i) {
        dev_staged_ids[base + i] = s_candidates[i].edge_id;
        dev_staged_dist2[base + i] = s_candidates[i].dist;
    }
}

//> Copy staged per-anchor rows into compact CSR arrays (no neighbor rediscovery).
__global__ void compact_staged_neighbors_kernel(
    int num_edges,
    int max_candidates,
    const int *dev_staged_counts,
    const int *dev_staged_ids,
    const float *dev_staged_dist2,
    const int *dev_neighbor_offsets,
    int *dev_neighbor_ids,
    float *dev_neighbor_dist2)
{
    const int te_idx = blockIdx.x;
    if (te_idx >= num_edges) {
        return;
    }

    const int count = dev_staged_counts[te_idx];
    const int base_out = dev_neighbor_offsets[te_idx];
    const int base_in = te_idx * max_candidates;
    for (int i = 0; i < count; ++i) {
        dev_neighbor_ids[base_out + i] = dev_staged_ids[base_in + i];
        dev_neighbor_dist2[base_out + i] = dev_staged_dist2[base_in + i];
    }
}

#endif // GPU_PREPROCESS_KERNELS_CUH
