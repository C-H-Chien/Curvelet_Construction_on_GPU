#ifndef GPU_DEV_EDGE_CHAIN_GROWTH_CUH
#define GPU_DEV_EDGE_CHAIN_GROWTH_CUH

#include <cmath>

#include "gpu_dev_utils.cuh"

__device__ __forceinline__ bool
slot_active_for_run(int f_run, float along)
{
    //> Style-2 direction filter: forward (f_run=0) keeps neighbors along +tangent; backward along -tangent.
    return (f_run == 0) ? (along > 0.f) : (along < 0.f);
}

__device__ __forceinline__ bool
check_curvelet_exist(
    int chain_sz,
    const unsigned *chain,
    int slots_per_anchor,
    int group_max_sz,
    const unsigned *edge_chain_target)
{
    const int row_w = group_max_sz + 1;
    for (int i = 0; i < slots_per_anchor; i++) {
        const unsigned stored_sz = edge_chain_target[i * row_w + group_max_sz];
        if (stored_sz == 0u) {
            return false;
        }
        if (static_cast<int>(stored_sz) != chain_sz) {
            continue;
        }
        bool match = true;
        for (int c = 0; c < chain_sz; c++) {
            if (chain[c] != edge_chain_target[i * row_w + c]) {
                match = false;
                break;
            }
        }
        if (match) {
            return true;
        }
    }
    return false;
}

//> Copy a pairwise (or other) bundle into the lane working grids (2 × bundle_cells).
__device__ __forceinline__ void
copy_bundle(
    int bundle_cells,
    float *dst_min_ks, float *dst_max_ks,
    const float *src_min_ks, const float *src_max_ks)
{
    for (int bidx = 0; bidx < bundle_cells; bidx++) {
        dst_min_ks[bidx] = src_min_ks[bidx];
        dst_max_ks[bidx] = src_max_ks[bidx];
    }
}

//> Interval intersection work ∩ cand per cell: new_min = max(mins), new_max = min(maxs).
//> Two-pass so a rejected candidate leaves the working bundle unchanged (no extra temp grids).
//> Returns true if the intersection is non-empty (any cell with new_max > new_min).
__device__ __forceinline__ bool
intersect_working_with_candidate(
    int bundle_cells,
    float *work_min_ks, float *work_max_ks,
    const float *cand_min_ks, const float *cand_max_ks)
{
    bool valid = false;
    for (int j = 0; j < bundle_cells; j++) {
        const float new_min = (work_min_ks[j] < cand_min_ks[j]) ? cand_min_ks[j] : work_min_ks[j];
        const float new_max = (work_max_ks[j] > cand_max_ks[j]) ? cand_max_ks[j] : work_max_ks[j];
        if (new_max > new_min) {
            valid = true;
            break;
        }
    }
    if (!valid) {
        return false;
    }
    for (int j = 0; j < bundle_cells; j++) {
        const float wmin = work_min_ks[j];
        const float wmax = work_max_ks[j];
        const float cmin = cand_min_ks[j];
        const float cmax = cand_max_ks[j];
        work_min_ks[j] = (wmin < cmin) ? cmin : wmin;
        work_max_ks[j] = (wmax > cmax) ? cmax : wmax;
    }
    return true;
}

__device__ __forceinline__ void
fill_curvelet_info_row(
    float *info_row,
    bool forward,
    float sx, float st,
    int curves_num_in_bundle_pixel, int curves_num_in_bundle_theta,
    const float *cmp_bundle_min_ks, const float *cmp_bundle_max_ks)
{
    //> Pick the valid (dx, dtheta) cell closest to the origin in the working bundle,
    //> then store only forward / k_max / k_min for that cell.
    float mind = 100.f;
    int mini = 0;
    int minj = 0;
    for (int ii = 0; ii < curves_num_in_bundle_pixel; ii++) {
        for (int jj = 0; jj < curves_num_in_bundle_theta; jj++) {
            const int bidx = ii * curves_num_in_bundle_theta + jj;
            if (cmp_bundle_max_ks[bidx] > cmp_bundle_min_ks[bidx]) {
                const float ddx = sx * (static_cast<float>(ii) - (static_cast<float>(curves_num_in_bundle_pixel) - 1.f) / 2.f);
                const float ddt = st * (static_cast<float>(jj) - (static_cast<float>(curves_num_in_bundle_theta) - 1.f) / 2.f);
                const float d = ddx * ddx + ddt * ddt;
                if (d < mind) {
                    mind = d;
                    mini = ii;
                    minj = jj;
                }
            }
        }
    }

    const int best_bidx = mini * curves_num_in_bundle_theta + minj;
    info_row[0] = forward ? 1.f : 0.f;
    info_row[1] = cmp_bundle_max_ks[best_bidx];
    info_row[2] = cmp_bundle_min_ks[best_bidx];

    // const float k_max = cmp_bundle_max_ks[best_bidx];
    // const float k_min = cmp_bundle_min_ks[best_bidx];
    // const float k = 0.5f * (k_max + k_min);
    // const float dx = sx * (static_cast<float>(mini) - (static_cast<float>(curves_num_in_bundle_pixel) - 1.f) / 2.f);
    // const float dt = st * (static_cast<float>(minj) - (static_cast<float>(curves_num_in_bundle_theta) - 1.f) / 2.f);
    // const float theta = To2Pi(ref_theta + dt);
    // const float pt_x = ref_pt_x - dx * sinf(theta);
    // const float pt_y = ref_pt_y + dx * cosf(theta);

    // float length = 0.f;
    // for (int i = 0; i + 1 < chain_sz; i++) {
    //     const int e1 = static_cast<int>(chain[i]);
    //     const int e2 = static_cast<int>(chain[i + 1]);
    //     const float x1 = dev_edges[e1 * sz_edge_data + 0];
    //     const float y1 = dev_edges[e1 * sz_edge_data + 1];
    //     const float x2 = dev_edges[e2 * sz_edge_data + 0];
    //     const float y2 = dev_edges[e2 * sz_edge_data + 1];
    //     const float ddx = x1 - x2;
    //     const float ddy = y1 - y2;
    //     length += sqrtf(ddx * ddx + ddy * ddy);
    // }

    // const float alpha3 = 1.f;
    // const float alpha4 = 1.f;
    // const float quality = (length > 0.f && chain_sz > 0)
    //     ? 2.f / (alpha3 * nrad / length + alpha4 * length / static_cast<float>(chain_sz))
    //     : 0.f;
    // info_row[1] = k_max;
    // info_row[2] = k_min;
    // info_row[3] = ref_theta;
    // info_row[4] = pt_x;
    // info_row[5] = pt_y;
    // info_row[6] = theta;
    // info_row[7] = k;
    // info_row[8] = length;
    // info_row[9] = quality;
}

//> Grow one seed chain into out_chain / out_len. Workspace is 2 × bundle_cells (working min + max).
//> On success with a grown chain, working grids hold the final bundle (for curvelet info).

//> Initialize a seed: direction/hyp filters, write anchor into chain, copy seed pairwise bundle into working.
//> Returns false if this seed is inactive for f_run.
__device__ __forceinline__ bool
grow_seed_init(
    int f_run,
    int seed_idx,
    int row_base,
    size_t anchor_bundle_base,
    int bundle_cells,
    int sz_edge_data,
    float anchor_edge_x,
    float anchor_edge_y,
    float anchor_cos,
    float anchor_sin,
    unsigned anchor_id,
    const float *dev_edges,
    const int *dev_neighbor_list,
    const float *dev_bundle_min_ks,
    const float *dev_bundle_max_ks,
    const unsigned char *dev_is_bundle_geometrically_valid,
    float *work_min_ks,
    float *work_max_ks,
    unsigned *out_chain,
    int *out_len)
{
    *out_len = 0;

    const int neighbor_id = dev_neighbor_list[row_base + seed_idx];
    if (neighbor_id < 0) {
        return false;
    }

    //> Directional filter: check if the seed neighbor is active for the current run (f_run == 0: forward, f_run == 1: backward)
    const size_t hyp_idx = static_cast<size_t>(row_base) + static_cast<size_t>(seed_idx);
    if (!dev_is_bundle_geometrically_valid[hyp_idx]) {
        return false;
    }
    const float *nbr = dev_edges + neighbor_id * sz_edge_data;
    const float dir = (nbr[0] - anchor_edge_x) * anchor_cos + (nbr[1] - anchor_edge_y) * anchor_sin;
    if (!slot_active_for_run(f_run, dir)) {
        return false;
    }

    out_chain[0] = anchor_id;
    *out_len = 1;

    //> Initialize working bundle from the seed pairwise curve bundle
    const float *seed_min = dev_bundle_min_ks + anchor_bundle_base + static_cast<size_t>(seed_idx) * static_cast<size_t>(bundle_cells);
    const float *seed_max = dev_bundle_max_ks + anchor_bundle_base + static_cast<size_t>(seed_idx) * static_cast<size_t>(bundle_cells);
    copy_bundle(bundle_cells, work_min_ks, work_max_ks, seed_min, seed_max);
    return true;
}

//> Consider one remaining neighbor for an active seed. cand_min/max are the pairwise grids for that slot.
//> Returns true if the chain hit group_max_sz and growth should stop.
__device__ __forceinline__ bool
grow_seed_consider_neighbor(
    int f_run,
    int seed_idx,
    int neighbor_remain_idx,
    int row_base,
    int bundle_cells,
    int group_max_sz,
    int sz_edge_data,
    float anchor_edge_x,
    float anchor_edge_y,
    float anchor_cos,
    float anchor_sin,
    const float *dev_edges,
    const int *dev_neighbor_list,
    const unsigned char *dev_is_bundle_geometrically_valid,
    const float *cand_min,
    const float *cand_max,
    float *work_min_ks,
    float *work_max_ks,
    unsigned *out_chain,
    int *edge_chain_nbr_idx)
{
    const int remain_id = dev_neighbor_list[row_base + neighbor_remain_idx];
    if (remain_id < 0) {
        return false;
    }

    const size_t hyp_remain = static_cast<size_t>(row_base) + static_cast<size_t>(neighbor_remain_idx);
    if (!dev_is_bundle_geometrically_valid[hyp_remain]) {
        return false;
    }

    //> Directional filter: check if the remaining neighbor is active for the current run (f_run == 0: forward, f_run == 1: backward)
    const float *nbr = dev_edges + remain_id * sz_edge_data;
    const float dir = (nbr[0] - anchor_edge_x) * anchor_cos + (nbr[1] - anchor_edge_y) * anchor_sin;
    if (!slot_active_for_run(f_run, dir)) {
        return false;
    }

    if (seed_idx == neighbor_remain_idx) {
        //> Style-2: ref_first == true -> append without re-intersecting the seed bundle
        if (*edge_chain_nbr_idx < group_max_sz) {
            out_chain[(*edge_chain_nbr_idx)++] = static_cast<unsigned>(remain_id);
        }
        return (*edge_chain_nbr_idx >= group_max_sz);
    }

    //> Stream cand against working; commit in place only if intersection is non-empty
    if (intersect_working_with_candidate(bundle_cells, work_min_ks, work_max_ks, cand_min, cand_max)) {
        if (*edge_chain_nbr_idx < group_max_sz) {
            out_chain[(*edge_chain_nbr_idx)++] = static_cast<unsigned>(remain_id);
        }
    }

    return (*edge_chain_nbr_idx >= group_max_sz);
}

__device__ __forceinline__ void
grow_seed_chain(
    int f_run,
    int seed_idx,
    int num_of_neighbors,
    int row_base,
    size_t anchor_bundle_base,
    int bundle_cells,
    int group_max_sz,
    int sz_edge_data,
    float anchor_edge_x,
    float anchor_edge_y,
    float anchor_cos,
    float anchor_sin,
    unsigned anchor_id,
    const float *dev_edges,
    const int *dev_neighbor_list,
    const float *dev_bundle_min_ks,
    const float *dev_bundle_max_ks,
    const unsigned char *dev_is_bundle_geometrically_valid,
    float *work_min_ks,
    float *work_max_ks,
    unsigned *out_chain,
    int *out_len)
{
    if (!grow_seed_init(
            f_run, seed_idx, row_base, anchor_bundle_base, bundle_cells, sz_edge_data,
            anchor_edge_x, anchor_edge_y, anchor_cos, anchor_sin, anchor_id,
            dev_edges, dev_neighbor_list,
            dev_bundle_min_ks, dev_bundle_max_ks, dev_is_bundle_geometrically_valid,
            work_min_ks, work_max_ks,
            out_chain, out_len)) {
        return;
    }

    //> Neighbors are already in ascending distance to the anchor
    for (int neighbor_remain_idx = 0; neighbor_remain_idx < num_of_neighbors; neighbor_remain_idx++) {
        const float *cand_min = dev_bundle_min_ks + anchor_bundle_base +
            static_cast<size_t>(neighbor_remain_idx) * static_cast<size_t>(bundle_cells);
        const float *cand_max = dev_bundle_max_ks + anchor_bundle_base +
            static_cast<size_t>(neighbor_remain_idx) * static_cast<size_t>(bundle_cells);

        if (grow_seed_consider_neighbor(
                f_run, seed_idx, neighbor_remain_idx, row_base,
                bundle_cells, group_max_sz, sz_edge_data,
                anchor_edge_x, anchor_edge_y, anchor_cos, anchor_sin,
                dev_edges, dev_neighbor_list, dev_is_bundle_geometrically_valid,
                cand_min, cand_max,
                work_min_ks, work_max_ks,
                out_chain, out_len)) {
            break;
        }
    }
}

__device__ __forceinline__ void
record_accepted_chain(
    unsigned anchor_id,
    int f_run,
    int chain_sz,
    const unsigned *chain,
    int group_max_sz,
    int chain_width,
    int max_per_anchor,
    float sx, float st,
    int curves_num_in_bundle_pixel, int curves_num_in_bundle_theta,
    const float *working_min_ks, const float *working_max_ks,
    unsigned *dev_edge_chain_final,
    unsigned *anchor_count,
    float *dev_curvelet_info)
{
    if (*anchor_count >= static_cast<unsigned>(max_per_anchor)) {
        return;
    }

    const unsigned row = static_cast<unsigned>(anchor_id) * static_cast<unsigned>(max_per_anchor) + (*anchor_count);
    (*anchor_count)++;

    unsigned *out_row = dev_edge_chain_final + static_cast<size_t>(row) * static_cast<size_t>(chain_width);
    for (int j = 0; j < chain_width; j++) {
        out_row[j] = 0;
    }
    out_row[0] = anchor_id + 1u;
    for (int i = 0; i < chain_sz && (i + 1) < chain_width; i++) {
        out_row[i + 1] = chain[i] + 1u;
    }

    float *info_row = dev_curvelet_info + static_cast<size_t>(row) * 10u;
    fill_curvelet_info_row(
        info_row, f_run == 0,
        sx, st,
        curves_num_in_bundle_pixel, curves_num_in_bundle_theta,
        working_min_ks, working_max_ks);

    (void)group_max_sz;
}

#endif // GPU_DEV_EDGE_CHAIN_GROWTH_CUH
