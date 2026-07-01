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
    SpatialIndexDevice spatial,
    float rad_sqr,
    unsigned neighbor_radius,
    int *dev_neighbor_counts,
    unsigned *dev_max_num_of_neighbors)
{
    const int anchor_edge_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (anchor_edge_idx >= num_edges) {
        return;
    }

    //> Do a binary search to identify valid neighbor edges. Here we only count the amount of neighbor edges.
    const int count = discover_neighbors_by_binary_search(anchor_edge_idx, spatial.dev_edges, spatial, rad_sqr, neighbor_radius, nullptr, 0, true);

    dev_neighbor_counts[anchor_edge_idx] = count;
    if (count > 0) {
        atomicMax(dev_max_num_of_neighbors, static_cast<unsigned>(count));
    }
}

__global__ void fill_csr_neighbors_fixed_kernel(
    int num_edges,
    unsigned max_num_of_neighbors,
    SpatialIndexDevice spatial,
    float rad_sqr,
    unsigned neighbor_radius,
    const int *dev_neighbor_counts,
    const int *dev_neighbor_offsets,
    int *dev_neighbor_ids,
    float *dev_neighbor_dist2)
{
    extern __shared__ NeighborCandidate s_candidates[];

    const int anchor_edge_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (anchor_edge_idx >= num_edges) {
        return;
    }

    //> Each thread owns max_num_of_neighbors slots in the dynamic shared memory.
    NeighborCandidate *my_candidates = s_candidates + threadIdx.x * max_num_of_neighbors;

    const int neighbor_count = dev_neighbor_counts[anchor_edge_idx];

    int count = discover_neighbors_by_binary_search(
        anchor_edge_idx, spatial.dev_edges, spatial, rad_sqr, neighbor_radius,
        my_candidates, neighbor_count);
    sort_neighbors_by_distance(my_candidates, count);

    const int base = dev_neighbor_offsets[anchor_edge_idx];
    for (int i = 0; i < count; ++i) {
        dev_neighbor_ids[base + i] = my_candidates[i].edge_id;
        dev_neighbor_dist2[base + i] = my_candidates[i].dist;
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

    const int anchor_edge_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (anchor_edge_idx >= num_edges) {
        return;
    }

    NeighborCandidate *my_candidates = s_candidates + threadIdx.x * max_candidates;

    const int raw_count = discover_neighbors_by_binary_search(
        anchor_edge_idx, spatial.dev_edges, spatial, rad_sqr, neighbor_radius,
        my_candidates, max_candidates);
    sort_neighbors_by_distance(my_candidates, raw_count);

    int count = raw_count;
    if (raw_count >= max_candidates) {
        count = max_candidates;
        atomicAdd(dev_truncated_anchors, 1u);
    }

    dev_staged_counts[anchor_edge_idx] = count;
    if (count > 0) {
        atomicMax(dev_max_num_of_neighbors, static_cast<unsigned>(count));
    }

    const int base = anchor_edge_idx * max_candidates;
    for (int i = 0; i < count; ++i) {
        dev_staged_ids[base + i] = my_candidates[i].edge_id;
        dev_staged_dist2[base + i] = my_candidates[i].dist;
    }
    for (int i = count; i < max_candidates; ++i) {
        dev_staged_ids[base + i] = -1;
        dev_staged_dist2[base + i] = 0.f;
    }
}

//> Fixed-row: one warp per anchor discovers 7x7 cells in parallel, sorts, writes neighbor_list[row + k].
__global__ void discover_fixed_row_warp_kernel(
    int num_edges,
    int max_candidates,
    int slots_per_anchor,
    int warps_per_block,
    SpatialIndexDevice spatial,
    float rad_sqr,
    unsigned neighbor_radius,
    int *dev_neighbor_counts,
    int *dev_neighbor_list,
    float *dev_neighbor_dist2_row,
    unsigned *dev_max_num_of_neighbors,
    unsigned *dev_truncated_anchors)
{
    extern __shared__ unsigned char raw_shmem[];

    const int lane = threadIdx.x & 31;
    const int warp_id = threadIdx.x >> 5;
    const int anchor_edge_idx = blockIdx.x * warps_per_block + warp_id;

    const size_t buf_bytes = static_cast<size_t>(warps_per_block) * static_cast<size_t>(max_candidates) * sizeof(NeighborCandidate);
    NeighborCandidate *warp_bufs = reinterpret_cast<NeighborCandidate *>(raw_shmem);
    int *warp_counts = reinterpret_cast<int *>(raw_shmem + buf_bytes);

    NeighborCandidate *my_buf = warp_bufs + static_cast<size_t>(warp_id) * static_cast<size_t>(max_candidates);

    if (anchor_edge_idx < num_edges) {
        if (lane == 0) {
            warp_counts[warp_id] = 0;
        }
        __syncwarp();

        const float te_x = spatial.dev_edges[anchor_edge_idx * kEdgeFields + 0];
        const float te_y = spatial.dev_edges[anchor_edge_idx * kEdgeFields + 1];
        const int cx = static_cast<int>(lrintf(te_x));
        const int cy = static_cast<int>(lrintf(te_y));
        const int side = 2 * static_cast<int>(neighbor_radius) + 1;
        const int num_cells = side * side;

        for (int cell = lane; cell < num_cells; cell += 32) {
            const int cell_dy = cell / side;
            const int cell_dx = cell % side;
            const int px = cx - static_cast<int>(neighbor_radius) + cell_dx;
            const int py = cy - static_cast<int>(neighbor_radius) + cell_dy;
            discover_cell_neighbors_warp_lane(
                anchor_edge_idx, px, py, te_x, te_y,
                spatial.dev_edges, spatial, rad_sqr,
                my_buf, &warp_counts[warp_id], max_candidates);
        }
        __syncwarp();

        int count = 0;
        if (lane == 0) {
            const int raw_count = warp_counts[warp_id];
            count = raw_count;
            if (raw_count >= max_candidates) {
                count = max_candidates;
                atomicAdd(dev_truncated_anchors, 1u);
            }
            sort_neighbors_by_distance(my_buf, count);
            dev_neighbor_counts[anchor_edge_idx] = count;
            if (count > 0) {
                atomicMax(dev_max_num_of_neighbors, static_cast<unsigned>(count));
            }
            warp_counts[warp_id] = count;
        }
        __syncwarp();
        count = warp_counts[warp_id];

        const int row = anchor_edge_idx * slots_per_anchor;
        if (lane < count) {
            dev_neighbor_list[row + lane] = my_buf[lane].edge_id;
            dev_neighbor_dist2_row[row + lane] = my_buf[lane].dist;
        }
        for (int k = count + lane; k < slots_per_anchor; k += 32) {
            dev_neighbor_list[row + k] = -1;
            dev_neighbor_dist2_row[row + k] = 0.f;
        }
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
    const int anchor_edge_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (anchor_edge_idx >= num_edges) {
        return;
    }

    const int count = dev_staged_counts[anchor_edge_idx];
    const int base_out = dev_neighbor_offsets[anchor_edge_idx];
    const int base_in = anchor_edge_idx * max_candidates;
    for (int i = 0; i < count; ++i) {
        dev_neighbor_ids[base_out + i] = dev_staged_ids[base_in + i];
        dev_neighbor_dist2[base_out + i] = dev_staged_dist2[base_in + i];
    }
}

#endif // GPU_PREPROCESS_KERNELS_CUH
