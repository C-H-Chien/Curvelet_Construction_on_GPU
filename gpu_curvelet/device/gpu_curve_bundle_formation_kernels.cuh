#ifndef GPU_CURVELET_KERNELS_CUH
#define GPU_CURVELET_KERNELS_CUH

#include "device/gpu_dev_common.cuh"
#include "gpu_dev_curve_bundle.cuh"

//> curve bundle formation
//> Configuration: one warp per anchor_id edge, one thread per anchor_id-neighbor edge pair.
//> Direction filtering (forward/backward) is deferred to chain growth; pairwise transport is direction-agnostic.
__global__ void curve_bundle_formation_kernel(
    int num_edges,
    int slots_per_anchor,
    int bundle_cells,
    int curves_num_in_bundle_pixel,
    int curves_num_in_bundle_theta,
    int warps_per_block,
    float dx,
    float dt,
    float sx,
    float st,
    float max_k,
    int sz_edge_data,
    const float *dev_edges,
    const int *dev_neighbor_list,
    const int *dev_neighbor_counts,
    float *dev_bundle_min_ks,
    float *dev_bundle_max_ks,
    unsigned char *dev_hyp_look_edge,
    unsigned *dev_valid_pair_count)
{
    extern __shared__ float warp_anchor_trig[];

    const int lane = threadIdx.x & 31;
    const int warp_id = threadIdx.x >> 5;

    //> Retrieve the anchor edge id corresponding to dev_edges
    const int anchor_id = blockIdx.x * warps_per_block + warp_id;

    //> Retrieve the anchor edge data from dev_edges
    const int anchor_offset_in_dev_edges = anchor_id * sz_edge_data;
    const float anchor_edge_x            = dev_edges[anchor_offset_in_dev_edges + 0];
    const float anchor_edge_y            = dev_edges[anchor_offset_in_dev_edges + 1];
    const float anchor_edge_theta        = dev_edges[anchor_offset_in_dev_edges + 2];

    //> One cos/sinf pair per warp in shared memory; lane 0 writes, all lanes read after sync.
    const int trig_base = warp_id * 2;
    if (lane == 0) {
        warp_anchor_trig[trig_base + 0] = cosf(anchor_edge_theta);
        warp_anchor_trig[trig_base + 1] = sinf(anchor_edge_theta);
    }
    __syncwarp();

    const int num_of_neighbor_edges = dev_neighbor_counts[anchor_id];
    const int row_base = anchor_id * slots_per_anchor;
    const size_t anchor_bundle_base = static_cast<size_t>(anchor_id) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);

    //> Loop over all neighbor edges
    for (int slot = lane; slot < slots_per_anchor; slot += 32) {
        const size_t hyp_idx = static_cast<size_t>(anchor_id) * static_cast<size_t>(slots_per_anchor) + static_cast<size_t>(slot);
        float *slot_min = dev_bundle_min_ks + anchor_bundle_base + static_cast<size_t>(slot) * static_cast<size_t>(bundle_cells);
        float *slot_max = dev_bundle_max_ks + anchor_bundle_base + static_cast<size_t>(slot) * static_cast<size_t>(bundle_cells);

        if (slot >= num_of_neighbor_edges) {
            dev_hyp_look_edge[hyp_idx] = 0;
            continue;
        }

        //> Retrieve the neighbor edge data
        const int neighbor_id                  = dev_neighbor_list[row_base + slot];
        if (neighbor_id < 0) {
            dev_hyp_look_edge[hyp_idx] = 0;
            continue;
        }
        const int neighbor_offset_in_dev_edges = neighbor_id * sz_edge_data;
        const float neighbor_edge_x            = dev_edges[neighbor_offset_in_dev_edges + 0];
        const float neighbor_edge_y            = dev_edges[neighbor_offset_in_dev_edges + 1];
        const float neighbor_edge_theta        = dev_edges[neighbor_offset_in_dev_edges + 2];

        //> Circular Arc Bundle Transport from Proposition 4
        const bool valid = get_circular_arc_bundle_transport(
            dx, dt, sx, st, max_k,
            curves_num_in_bundle_pixel, curves_num_in_bundle_theta,
            slot_min, slot_max,
            anchor_edge_x, anchor_edge_y,
            warp_anchor_trig[trig_base + 0], warp_anchor_trig[trig_base + 1],
            neighbor_edge_x, neighbor_edge_y, neighbor_edge_theta);

        dev_hyp_look_edge[hyp_idx] = valid ? 1 : 0;
        if (valid && dev_valid_pair_count != nullptr) {
            atomicAdd(dev_valid_pair_count, 1u);
        }
    }
}

#endif // GPU_CURVELET_KERNELS_CUH
