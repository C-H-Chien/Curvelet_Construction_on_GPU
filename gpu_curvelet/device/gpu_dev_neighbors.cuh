#ifndef GPU_DEV_NEIGHBORS_CUH
#define GPU_DEV_NEIGHBORS_CUH

#include "gpu_dev_common.cuh"

struct NeighborCandidate {
    float dist;
    int edge_id;
};

struct SpatialIndexDevice {
    float *dev_edges = nullptr;
    long long *dev_unique_cells = nullptr;
    int *dev_cell_starts = nullptr;
    int *dev_cell_counts = nullptr;
    int *dev_sorted_edge_ids = nullptr;
    int num_cells = 0;
};

__device__ int discover_neighbors_by_binary_search(
    int te_idx,
    const float *dev_edges,
    const SpatialIndexDevice &spatial,
    float rad_sqr,
    unsigned neighbor_radius,
    NeighborCandidate *out,
    int max_out)
{
    const float te_x = dev_edges[te_idx * kEdgeFields + 0];
    const float te_y = dev_edges[te_idx * kEdgeFields + 1];
    const int cx = static_cast<int>(lrintf(te_x));
    const int cy = static_cast<int>(lrintf(te_y));

    int count = 0;

    //> Loop over the neighbor radius in x-direction
    for (int py = cy - static_cast<int>(neighbor_radius); py <= cy + static_cast<int>(neighbor_radius); ++py) {

        //> Loop over the neighbor radius in y-direction
        for (int px = cx - static_cast<int>(neighbor_radius); px <= cx + static_cast<int>(neighbor_radius); ++px) {
            
            const long long key = pack_cell_key(px, py);
            const int cell_idx = lower_bound_cell_keys(
                spatial.dev_unique_cells, spatial.num_cells, key);
            if (cell_idx >= spatial.num_cells ||
                spatial.dev_unique_cells[cell_idx] != key) {
                continue;
            }

            //> Fetch the offset
            const int start = spatial.dev_cell_starts[cell_idx];
            const int end   = start + spatial.dev_cell_counts[cell_idx];

            //> Loop over all the edges in the current cell
            for (int k = start; k < end; ++k) {
                const int ne_idx = spatial.dev_sorted_edge_ids[k];
                if (ne_idx == te_idx) {
                    continue;
                }

                //> Fetch the neighbor edge coordinates
                const float ne_x = dev_edges[ne_idx * kEdgeFields + 0];
                const float ne_y = dev_edges[ne_idx * kEdgeFields + 1];
                const float dist = sq_dist_dev(te_x, te_y, ne_x, ne_y);
                if (dist > rad_sqr) {
                    continue;
                }

                if (count < max_out) {
                    out[count].dist = dist;
                    out[count].edge_id = ne_idx;
                    ++count;
                }
            }
        }
    }
    return count;
}

__device__ void sort_neighbors_by_distance(NeighborCandidate *candidates, int count)
{
    for (int i = 1; i < count; ++i) {
        NeighborCandidate key = candidates[i];
        int j = i - 1;

        //> Neighbors are sorted by distance in an ascending order
        while (j >= 0 && candidates[j].dist > key.dist) {
            candidates[j + 1] = candidates[j];
            --j;
        }
        candidates[j + 1] = key;
    }
}

#endif // GPU_DEV_NEIGHBORS_CUH
