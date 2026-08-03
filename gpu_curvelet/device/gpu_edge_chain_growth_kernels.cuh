#ifndef GPU_EDGE_CHAIN_GROWTH_KERNELS_CUH
#define GPU_EDGE_CHAIN_GROWTH_KERNELS_CUH

#include <cmath>

#include "gpu_dev_edge_chain_growth.cuh"
#include "gpu_dev_utils.cuh"

//> Shared-memory cache modes for warp growth (host selects based on size limits).
//> 0 = none (global bundles + global lane workspaces)
//> 1 = lane workspaces in shared; pairwise bundles stay in global
enum : int {
    kWarpSmemNone = 0,
    kWarpSmemLaneWs = 1
};

//> One warp per anchor: Phase 1
//> Lanes grow different seeds in parallel for both forward and backward directions,
//> writing into separate per-direction candidate / working-bundle buffers.
//>
//> Dynamic shared memory (per warp in the block), when enabled:
//>   [32 lanes × 2*bundle_cells growth workspace] // mode 1 (working min+max)
//>
//> Float scratch per anchor (global; durable across phase 1→2):
//>   [2 × slots × bundle_cells seed working min]   // f_run = 0,1
//>   [2 × slots × bundle_cells seed working max]   // f_run = 0,1
//>   [32 lanes × 2*bundle_cells growth workspace]  // only when smem_mode == 0
//> UInt scratch per anchor:
//>   [2 × slots × chain_width candidate chains]    // f_run = 0,1; length in last column
//>   [slots × chain_width dedup table]             // used by phase-2 kernel
__global__ void grow_edge_chains_warp_kernel(
    int num_edges,
    int slots_per_anchor,
    int bundle_cells,
    int group_max_sz,
    int chain_width,
    int sz_edge_data,
    int warps_per_block,
    int smem_mode,
    size_t scratch_floats_per_anchor,
    size_t scratch_uints_per_anchor,
    const float *dev_edges,
    const int *dev_neighbor_list,
    const int *dev_neighbor_counts,
    const float *dev_bundle_min_ks,
    const float *dev_bundle_max_ks,
    const unsigned char *dev_is_bundle_geometrically_valid,
    float *dev_scratch_f,
    unsigned *dev_scratch_u)
{
    extern __shared__ float smem[];

    const int lane = threadIdx.x & 31;
    const int warp_id = threadIdx.x >> 5;
    const int anchor_id = static_cast<int>(blockIdx.x * warps_per_block + warp_id);
    if (anchor_id >= num_edges) {
        return;
    }

    const int num_of_neighbors = dev_neighbor_counts[anchor_id];
    const int row_base = anchor_id * slots_per_anchor;
    const size_t anchor_bundle_base = static_cast<size_t>(anchor_id) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);

    //> Retrieve anchor edge information
    const float *anchor_edge = dev_edges + anchor_id * sz_edge_data;
    const float anchor_edge_x = anchor_edge[0];
    const float anchor_edge_y = anchor_edge[1];
    const float anchor_orient = anchor_edge[2];
    const float anchor_cos = cosf(anchor_orient);
    const float anchor_sin = sinf(anchor_orient);

    const size_t lane_workspace_floats = static_cast<size_t>(2) * static_cast<size_t>(bundle_cells);
    const bool lane_ws_in_smem = (smem_mode == kWarpSmemLaneWs);
    const size_t floats_per_warp = lane_ws_in_smem ? (static_cast<size_t>(32) * lane_workspace_floats)
        : 0;
    float *warp_smem = (floats_per_warp > 0) ? (smem + static_cast<size_t>(warp_id) * floats_per_warp) : nullptr;

    //> candidate_chains[f_run]: slots * chain_width; length at column group_max_sz
    unsigned *candidate_chains = dev_scratch_u + static_cast<size_t>(anchor_id) * scratch_uints_per_anchor;

    //> float stores curvature grids (bundle min/max while comparing/intersecting)
    //> unsigned stores edge IDs in chains
    float *seed_working_min = dev_scratch_f + static_cast<size_t>(anchor_id) * scratch_floats_per_anchor;
    float *seed_working_max = seed_working_min + static_cast<size_t>(2) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);

    float *lane_ws_base = nullptr;
    if (lane_ws_in_smem) {
        lane_ws_base = warp_smem;
    } 
    else {
        //> Fallback: lane workspaces live after seed-working grids in global scratch
        lane_ws_base = seed_working_max + static_cast<size_t>(2) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);
    }

    //> if shared memory is used, the following pointers are pointing to the shared memory workspace
    float *lane_ws = lane_ws_base + static_cast<size_t>(lane) * lane_workspace_floats;
    float *work_min_ks = lane_ws;
    float *work_max_ks = lane_ws + bundle_cells;
    
    //> f_run = 0: forward growth
    //> f_run = 1: backward growth
    // #pragma unroll 2
    for (int f_run = 0; f_run < 2; f_run++) {
        //> Each lane grows a disjoint set of seeds into candidate_chains[f_run]
        unsigned *candidate_chains_frun = candidate_chains + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(chain_width);
        float *seed_working_min_frun = seed_working_min + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);
        float *seed_working_max_frun = seed_working_max + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);

        //> Loop over all seeds in the neighbor slots of the anchor
        for (int seed_idx = lane; seed_idx < slots_per_anchor; seed_idx += 32) {
            unsigned *cand_row = candidate_chains_frun + static_cast<size_t>(seed_idx) * static_cast<size_t>(chain_width);
            for (int j = 0; j < chain_width; j++) {
                cand_row[j] = 0;
            }

            float *out_min = seed_working_min_frun + static_cast<size_t>(seed_idx) * static_cast<size_t>(bundle_cells);
            float *out_max = seed_working_max_frun + static_cast<size_t>(seed_idx) * static_cast<size_t>(bundle_cells);

            if (seed_idx >= num_of_neighbors) {
                continue;
            }

            //> Grow directly into cand_row[0...group_max_sz).
            int chain_len = 0;
            grow_seed_chain(
                f_run, seed_idx, num_of_neighbors, row_base, anchor_bundle_base,
                bundle_cells, group_max_sz, sz_edge_data,
                anchor_edge_x, anchor_edge_y, anchor_cos, anchor_sin,
                static_cast<unsigned>(anchor_id),
                dev_edges, dev_neighbor_list,
                dev_bundle_min_ks, dev_bundle_max_ks, dev_is_bundle_geometrically_valid,
                work_min_ks, work_max_ks,
                cand_row, &chain_len);

            //> The curvelet length is stored at cand_row[group_max_sz]
            cand_row[group_max_sz] = static_cast<unsigned>(chain_len);

            if (chain_len > 0) {
                for (int b = 0; b < bundle_cells; b++) {
                    out_min[b] = work_min_ks[b];
                    out_max[b] = work_max_ks[b];
                }
            }
        }
        __syncwarp();
    }

    (void)scratch_uints_per_anchor;
    (void)anchor_orient;
}

//> One thread per anchor: Phase 2
//> Reads both-direction candidate chains + working bundles from warp growth scratch,
//> then serializes exact-match dedup and output recording (f_run = 0 then f_run = 1).
__global__ void dedup_record_edge_chains_kernel(
    int num_edges,
    int slots_per_anchor,
    int bundle_cells,
    int curves_num_in_bundle_pixel,
    int curves_num_in_bundle_theta,
    int group_max_sz,
    int max_per_anchor,
    int chain_width,
    float sx,
    float st,
    size_t scratch_floats_per_anchor,
    size_t scratch_uints_per_anchor,
    const int *dev_neighbor_counts,
    float *dev_scratch_f,
    unsigned *dev_scratch_u,
    unsigned *dev_edge_chain_final,
    unsigned *dev_anchor_chain_count,
    float *dev_curvelet_info)
{
    const int anchor_id = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (anchor_id >= num_edges) {
        return;
    }

    const int num_of_neighbors = dev_neighbor_counts[anchor_id];

    float *scratch_f = dev_scratch_f + static_cast<size_t>(anchor_id) * scratch_floats_per_anchor;
    unsigned *scratch_u = dev_scratch_u + static_cast<size_t>(anchor_id) * scratch_uints_per_anchor;

    //> Must match grow_edge_chains_warp_kernel scratch layout:
    //> [2 × slots × cells seed working min][2 × slots × cells seed working max][(optional lane WS)]
    float *seed_working_min = scratch_f;
    float *seed_working_max = seed_working_min + static_cast<size_t>(2) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);

    //> candidate_chains for both directions, then dedup table
    unsigned *candidate_chains = scratch_u;
    unsigned *edge_chain_target = scratch_u + static_cast<size_t>(2) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(chain_width);

    unsigned &anchor_count = dev_anchor_chain_count[anchor_id];
    anchor_count = 0;

    const int target_row_w = chain_width;

    //> f_run = 0: forward; f_run = 1: backward; dedup tables are independent per direction
    // #pragma unroll 2
    for (int f_run = 0; f_run < 2; f_run++) {
        //> Serialize exact-match dedup and output recording (seed order)
        int edge_chain_nbr_idx = 0;
        for (int i = 0; i < slots_per_anchor * target_row_w; i++) {
            edge_chain_target[i] = 0;
        }

        //> Get the candidate chains for the current growth direction
        unsigned *candidate_chains_frun = candidate_chains + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(chain_width);
        //> Get the seed working min for the current growth direction
        float *seed_working_min_frun = seed_working_min + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);
        float *seed_working_max_frun = seed_working_max + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);

        //> Walks seeds in order
        for (int seed_idx = 0; seed_idx < num_of_neighbors; seed_idx++) {
            
            //> Get the candidate chain for the current seed
            unsigned *cand_row = candidate_chains_frun + static_cast<size_t>(seed_idx) * static_cast<size_t>(chain_width);
            const int chain_len = static_cast<int>(cand_row[group_max_sz]);

            //> Skip the chain if the length is less than 2
            if (chain_len <= 2) {
                continue;
            }

            //> If the curvelet does not exist in the dedup table, accept it
            if (check_curvelet_exist(chain_len, cand_row, slots_per_anchor, group_max_sz, edge_chain_target)) {
                continue;
            }

            if (edge_chain_nbr_idx < slots_per_anchor) {
                for (int c = 0; c < chain_len; c++) {
                    edge_chain_target[edge_chain_nbr_idx * target_row_w + c] = cand_row[c];
                }
                edge_chain_target[edge_chain_nbr_idx * target_row_w + group_max_sz] = static_cast<unsigned>(chain_len);
                edge_chain_nbr_idx++;
            }

            const float *work_min = seed_working_min_frun + static_cast<size_t>(seed_idx) * static_cast<size_t>(bundle_cells);
            const float *work_max = seed_working_max_frun + static_cast<size_t>(seed_idx) * static_cast<size_t>(bundle_cells);
            record_accepted_chain(
                static_cast<unsigned>(anchor_id), f_run, chain_len, cand_row,
                group_max_sz, chain_width, max_per_anchor,
                sx, st,
                curves_num_in_bundle_pixel, curves_num_in_bundle_theta,
                work_min, work_max,
                dev_edge_chain_final, &anchor_count, dev_curvelet_info);
        }
    }

    (void)scratch_uints_per_anchor;
}

#endif // GPU_EDGE_CHAIN_GROWTH_KERNELS_CUH
