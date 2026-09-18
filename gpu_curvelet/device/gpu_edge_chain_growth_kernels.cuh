#ifndef GPU_EDGE_CHAIN_GROWTH_KERNELS_CUH
#define GPU_EDGE_CHAIN_GROWTH_KERNELS_CUH

#include <cmath>

#include "gpu_dev_edge_chain_growth.cuh"
#include "gpu_dev_utils.cuh"

//> Shared-memory cache modes for warp growth (host selects based on size limits).
//> 0 = none (global bundles + global lane workspaces)
//> 1 = lane workspaces in shared; pairwise bundles stay in global
//> 2 = lane workspaces + per-anchor pairwise bundles in shared
//> 3 = tile: one cooperatively loaded pairwise tile + a small working-grid pool
enum : int {
    kWarpSmemNone = 0,
    kWarpSmemLaneWs = 1,
    kWarpSmemLaneWsAndBundles = 2,
    kWarpSmemTile = 3
};


//> One warp per anchor: Phase 1
//> Lanes grow different seeds in parallel for both forward and backward directions,
//> writing into separate per-direction candidate / working-bundle buffers.
//>
//> Dynamic shared memory (per warp in the block), when enabled:
//>   mode 1: [32 lanes × 2*bundle_cells growth workspace]
//>   mode 2: [32 lanes × 2*bundle_cells growth workspace]
//>           [slots × bundle_cells min][slots × bundle_cells max]
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
    const size_t bundle_grid_floats = static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);
    const bool lane_ws_in_smem = (smem_mode == kWarpSmemLaneWs || smem_mode == kWarpSmemLaneWsAndBundles);
    const bool bundles_in_smem = (smem_mode == kWarpSmemLaneWsAndBundles);

    size_t floats_per_warp = 0;
    if (smem_mode == kWarpSmemLaneWs) {
        floats_per_warp = static_cast<size_t>(32) * lane_workspace_floats;
    } 
    else if (smem_mode == kWarpSmemLaneWsAndBundles) {
        floats_per_warp = (bundle_grid_floats * 2) + (static_cast<size_t>(32) * lane_workspace_floats);
    }
    float *warp_smem = (floats_per_warp > 0) ? (smem + static_cast<size_t>(warp_id) * floats_per_warp) : nullptr;

    //> candidate_chains[f_run]: slots * chain_width; length at column group_max_sz
    unsigned *candidate_chains = dev_scratch_u + static_cast<size_t>(anchor_id) * scratch_uints_per_anchor;

    //> float stores curvature grids (bundle min/max while comparing/intersecting)
    //> unsigned stores edge IDs in chains
    float *seed_working_min = dev_scratch_f + static_cast<size_t>(anchor_id) * scratch_floats_per_anchor;
    float *seed_working_max = seed_working_min + static_cast<size_t>(2) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);

    //> Pairwise bundles for this anchor: shared cache or global base
    const float *anchor_bundle_min = nullptr;
    const float *anchor_bundle_max = nullptr;
    float *lane_ws_base = nullptr;

    if (bundles_in_smem) {
        float *smem_bundle_min = warp_smem;
        float *smem_bundle_max = smem_bundle_min + bundle_grid_floats;
        lane_ws_base = smem_bundle_max + bundle_grid_floats;

        //> Cooperative load of this anchor's pairwise min/max grids into shared memory
        const size_t n_cells = static_cast<size_t>(num_of_neighbors) * static_cast<size_t>(bundle_cells);
        for (size_t i = static_cast<size_t>(lane); i < n_cells; i += 32) {
            smem_bundle_min[i] = dev_bundle_min_ks[anchor_bundle_base + i];
            smem_bundle_max[i] = dev_bundle_max_ks[anchor_bundle_base + i];
        }
        __syncwarp();

        anchor_bundle_min = smem_bundle_min;
        anchor_bundle_max = smem_bundle_max;
    } 
    else {
        anchor_bundle_min = dev_bundle_min_ks + anchor_bundle_base;
        anchor_bundle_max = dev_bundle_max_ks + anchor_bundle_base;
        if (lane_ws_in_smem) {
            lane_ws_base = warp_smem;
        } 
        else {
            //> Fallback: lane workspaces live after seed-working grids in global scratch
            lane_ws_base = seed_working_max + static_cast<size_t>(2) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);
        }
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
                f_run, seed_idx, num_of_neighbors, row_base,
                bundle_cells, group_max_sz, sz_edge_data,
                anchor_edge_x, anchor_edge_y, anchor_cos, anchor_sin,
                static_cast<unsigned>(anchor_id),
                dev_edges, dev_neighbor_list,
                anchor_bundle_min, anchor_bundle_max, dev_is_bundle_geometrically_valid,
                work_min_ks, work_max_ks,
                cand_row, &chain_len);

            //> The curvelet length is stored at cand_row[group_max_sz]
            cand_row[group_max_sz] = static_cast<unsigned>(chain_len);

            //> Phase 2 only needs working grids for chains that pass the length gate
            if (chain_len > 2) {
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

//> One warp per anchor, streaming pairwise tiles: Phase 1 (--chain-smem-mode=tile)
//>
//> Inverts the seed loop so the warp loads each pairwise min/max grid once per wave
//> into a single shared tile, instead of caching every slot or re-reading it per seed.
//> Dynamic shared memory (per warp):
//>   [W × 2*bundle_cells working grids]
//>   [2*bundle_cells pairwise tile]
//>   [slots active-slot ids]   // after the float region, all warps
//> W default 8. If more than W neighbors pass hyp+direction, extra waves reload the tile.
//> Chain order matches style-2: compact active slots in ascending slot (distance) order,
//> append the seed without re-intersecting, stop a seed at group_max_sz.
//> Persists compact scratch: length column, and k_max/k_min if len > 2 (scratch-traffic cut).
__global__ void grow_edge_chains_tile_kernel(
    int num_edges,
    int slots_per_anchor,
    int bundle_cells,
    int group_max_sz,
    int chain_width,
    int sz_edge_data,
    int curves_num_in_bundle_pixel,
    int curves_num_in_bundle_theta,
    float sx,
    float st,
    int warps_per_block,
    int tile_workspaces,
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

    const float *anchor_edge = dev_edges + anchor_id * sz_edge_data;
    const float anchor_edge_x = anchor_edge[0];
    const float anchor_edge_y = anchor_edge[1];
    const float anchor_orient = anchor_edge[2];
    const float anchor_cos = cosf(anchor_orient);
    const float anchor_sin = sinf(anchor_orient);

    const size_t floats_per_warp = static_cast<size_t>(tile_workspaces + 1) * static_cast<size_t>(2) * static_cast<size_t>(bundle_cells);
    float *s_warp_f    = smem + static_cast<size_t>(warp_id) * floats_per_warp;
    float *s_tile_min  = s_warp_f + static_cast<size_t>(tile_workspaces) * static_cast<size_t>(2) * static_cast<size_t>(bundle_cells);
    float *s_tile_max  = s_tile_min + bundle_cells;
    int *active_slots  = reinterpret_cast<int *>(smem + static_cast<size_t>(warps_per_block) * floats_per_warp)
                        + static_cast<size_t>(warp_id) * static_cast<size_t>(slots_per_anchor);

    const float *src_min = dev_bundle_min_ks + anchor_bundle_base;
    const float *src_max = dev_bundle_max_ks + anchor_bundle_base;

    unsigned *candidate_chains = dev_scratch_u + static_cast<size_t>(anchor_id) * scratch_uints_per_anchor;
    float *seed_ks = dev_scratch_f + static_cast<size_t>(anchor_id) * scratch_floats_per_anchor;

    for (int f_run = 0; f_run < 2; f_run++) {
        unsigned *candidate_chains_frun = candidate_chains + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(chain_width);
        float *seed_ks_frun = seed_ks + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * 2u;

        for (int seed_idx = lane; seed_idx < num_of_neighbors; seed_idx += 32) {
            unsigned *cand_row = candidate_chains_frun + static_cast<size_t>(seed_idx) * static_cast<size_t>(chain_width);
            cand_row[group_max_sz] = 0u;
        }
        __syncwarp();

        //> First apply the directional filter so the "active slots" are neighbors that pass the filter
        int n_local = 0;
        if (lane == 0) {
            for (int slot = 0; slot < num_of_neighbors; slot++) {
                if (neighbor_slot_passes_grow_filter(
                        f_run, slot, row_base, sz_edge_data,
                        anchor_edge_x, anchor_edge_y, anchor_cos, anchor_sin,
                        dev_edges, dev_neighbor_list, dev_is_bundle_geometrically_valid)) {
                    active_slots[n_local++] = slot;
                }
            }
        }
        __syncwarp();
        //> broadcast the number of active slots to all lanes
        const int n_active = __shfl_sync(0xffffffff, n_local, 0);

        //> Loop over all active slots in waves
        const int ws_stride = 2 * bundle_cells;
        for (int wave_base = 0; wave_base < n_active; wave_base += tile_workspaces) {
            const int wave_n = (n_active - wave_base < tile_workspaces) ? (n_active - wave_base) : (tile_workspaces);

            int chain_len = 0;
            int my_seed_slot = -1;
            unsigned *cand_row = nullptr;
            float *s_work_min = nullptr;
            float *s_work_max = nullptr;
            bool seed_open = (lane < wave_n);

            for (int s = 0; s < wave_n; s++) {
                const int seed_slot = active_slots[wave_base + s];
                //> Cooperative load the pairwise min/max grids from global memory into the tile [s_tile_min, s_tile_max]
                coop_load_pairwise_bundle(seed_slot, bundle_cells, lane, src_min, src_max, s_tile_min, s_tile_max);
                __syncwarp();

                //> copy the tile [s_tile_min, s_tile_max] to the working grids [s_work_min, s_work_max]
                if (lane == s) {
                    my_seed_slot = seed_slot;
                    s_work_min = s_warp_f + static_cast<size_t>(lane) * static_cast<size_t>(ws_stride);
                    s_work_max = s_work_min + bundle_cells;
                    copy_bundle(bundle_cells, s_work_min, s_work_max, s_tile_min, s_tile_max);
                    cand_row = candidate_chains_frun + static_cast<size_t>(seed_slot) * static_cast<size_t>(chain_width);
                    cand_row[0] = static_cast<unsigned>(anchor_id);
                    chain_len = 1;
                }
                __syncwarp();
            }

            for (int a = 0; a < n_active; a++) {
                const int remain_slot = active_slots[a];
                //> Cooperative load the pairwise min/max grids from global memory again into the tile [s_tile_min, s_tile_max] for the remain (staging) slot
                coop_load_pairwise_bundle(remain_slot, bundle_cells, lane, src_min, src_max, s_tile_min, s_tile_max);
                __syncwarp();
                if (seed_open) {
                    const int remain_id = dev_neighbor_list[row_base + remain_slot];

                    //> If the remain slot is the same as the my_seed_slot, append the remain_id to the candidate chain
                    //> otherwise, intersect the working grids with the tile and append the remain_id to the candidate chain if the edge chain size permits
                    if (remain_slot == my_seed_slot) {
                        if (chain_len < group_max_sz) {
                            cand_row[chain_len++] = static_cast<unsigned>(remain_id);
                        }
                    }
                    else if (intersect_working_with_candidate(bundle_cells, s_work_min, s_work_max, s_tile_min, s_tile_max)) {
                        if (chain_len < group_max_sz) {
                            cand_row[chain_len++] = static_cast<unsigned>(remain_id);
                        }
                    }
                    if (chain_len >= group_max_sz) {
                        seed_open = false;
                    }
                }
                __syncwarp();
            }

            if (lane < wave_n) {
                cand_row[group_max_sz] = static_cast<unsigned>(chain_len);
                if (chain_len > 2) {
                    float k_max = 0.f;
                    float k_min = 0.f;
                    select_best_bundle_ks(
                        sx, st, curves_num_in_bundle_pixel, curves_num_in_bundle_theta,
                        work_min, work_max, k_max, k_min);
                    float *ks = seed_ks_frun + static_cast<size_t>(my_seed_slot) * 2u;
                    //> Write the active k_max and k_min to global memory
                    ks[0] = k_max;
                    ks[1] = k_min;
                }
            }
            __syncwarp();
        }
    }
}


//> One thread per anchor: Phase 2
//> Reads both-direction candidate chains from phase-1 scratch, then serializes
//> exact-match dedup and output recording (f_run = 0 then f_run = 1).
//> compact_ks=0: warp path stores full working min/max grids.
//> compact_ks=1: tile path stores only k_max/k_min per seed.
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
    int compact_ks,
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

    //> Must match phase-1 scratch layout:
    //>   compact_ks=0: [2 × slots × cells seed working min][2 × slots × cells seed working max][(optional lane WS)]
    //>   compact_ks=1: [2 × slots × 2 k_max/k_min]
    float *seed_working_min = scratch_f;
    float *seed_working_max = nullptr;
    float *seed_ks = scratch_f;
    if (!compact_ks) {
        seed_working_max = seed_working_min
                         + static_cast<size_t>(2) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);
    }

    //> candidate_chains for both directions, then dedup table
    unsigned *candidate_chains = scratch_u;
    unsigned *edge_chain_target = scratch_u + static_cast<size_t>(2) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(chain_width);

    unsigned &anchor_count = dev_anchor_chain_count[anchor_id];
    anchor_count = 0;

    const int target_row_w = chain_width;

    //> f_run = 0: forward; f_run = 1: backward; dedup tables are independent per direction
    for (int f_run = 0; f_run < 2; f_run++) {
        //> Serialize exact-match dedup and output recording (seed order)
        int edge_chain_nbr_idx = 0;
        if (compact_ks) {
            for (int i = 0; i < slots_per_anchor; i++) {
                edge_chain_target[i * target_row_w + group_max_sz] = 0;
            }
        }
        else {
            for (int i = 0; i < slots_per_anchor * target_row_w; i++) {
                edge_chain_target[i] = 0;
            }
        }

        //> Get the candidate chains for the current growth direction
        unsigned *candidate_chains_frun = candidate_chains + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(chain_width);
        float *seed_working_min_frun = nullptr;
        float *seed_working_max_frun = nullptr;
        float *seed_ks_frun = nullptr;
        if (compact_ks) {
            seed_ks_frun = seed_ks + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * 2u;
        }
        else {
            //> Get the seed working min/max for the current growth direction
            seed_working_min_frun = seed_working_min + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);
            seed_working_max_frun = seed_working_max + static_cast<size_t>(f_run) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);
        }

        //> Walks seeds in order
        for (int seed_idx = 0; seed_idx < num_of_neighbors; seed_idx++) {
            //> Get the candidate chain for the current seed
            unsigned *cand_row = candidate_chains_frun + static_cast<size_t>(seed_idx) * static_cast<size_t>(chain_width);

            //> Skip the chain if the length is less than 2
            const int chain_len = static_cast<int>(cand_row[group_max_sz]);
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

            if (compact_ks) {
                const float *ks = seed_ks_frun + static_cast<size_t>(seed_idx) * 2u;
                record_accepted_chain(
                    static_cast<unsigned>(anchor_id), f_run, chain_len, cand_row,
                    group_max_sz, chain_width, max_per_anchor,
                    ks[0], ks[1],
                    dev_edge_chain_final, &anchor_count, dev_curvelet_info);
            }
            else {
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
    }

    (void)scratch_uints_per_anchor;
}

#endif // GPU_EDGE_CHAIN_GROWTH_KERNELS_CUH
