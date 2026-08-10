#ifndef GPU_PREPROCESS_HPP
#define GPU_PREPROCESS_HPP

#include <cstddef>
#include <string>

#include "param_settings.hpp"

class CategoryProfiler;

// GPU preprocessing for the neighbor-search stage of curvelet construction.
//
// NeighborCSR (compressed sparse rows format):
//   Global edge attributes live once in dev_edges (num_edges x sz_edge_data).
//   Per-anchor neighbors are stored in CSR form:
//     neighbor_offsets[anchor] .. neighbor_offsets[anchor+1]  -> slice in neighbor_ids[]
//     neighbor_ids[k] is sorted by ascending squared distance to its anchor.
//     neighbor_dist2[k] holds that squared distance (optional but useful).
//
// Internally the spatial cell index is also CSR over pixel cells.

//> Neighbor graph after preprocessing (CSR or fixed-row layout).
struct GPUNeighborGraph {
    std::string layout = "fixed-row";

    int num_edges                 = 0;   //> number of edges
    int total_neighbor_pairs      = 0;   //> CSR: neighbor_offsets[num_edges]; sum of neighbor counts
    unsigned max_num_of_neighbors = 0;   //> max number of neighbors per anchor

    //> Global edge table on device (row-major: x, y, orientation, strength)
    float *dev_edges = nullptr;

    //> CSR layout: neighbors of anchor a in
    //>   neighbor_ids[ neighbor_offsets[a] .. neighbor_offsets[a+1] )
    int *dev_neighbor_offsets   = nullptr;    //> size num_edges + 1
    int *dev_neighbor_ids       = nullptr;    //> size total_neighbor_pairs
    float *dev_neighbor_dist2   = nullptr;    //> size total_neighbor_pairs

    //> Fixed-row layout (edgeLookList-style): row a starts at anchor * neighbor_slots_per_anchor
    int neighbor_slots_per_anchor = 0;        //> typically max_candidates
    int *dev_neighbor_list        = nullptr;  //> size num_edges * neighbor_slots_per_anchor, -1 padding
    float *dev_neighbor_dist2_row = nullptr;  //> same shape; unused slots are 0
    int *dev_neighbor_counts      = nullptr;  //> size num_edges; valid neighbors per anchor

    //> (Optional) host mirrors
    float *host_edges              = nullptr;
    int *host_neighbor_offsets     = nullptr;
    int *host_neighbor_ids         = nullptr;
    float *host_neighbor_dist2     = nullptr;
    int *host_neighbor_list        = nullptr;
    float *host_neighbor_dist2_row = nullptr;
    int *host_neighbor_counts      = nullptr;
};

//> Build neighbor data on the GPU from a host-side third-order edge array.
//> Optional profiler aggregates timings across pipeline phases (spatial index, neighbor graph, etc.).
bool gpu_preprocess_build(
    const CurveletParams &params,
    int device_id,
    int num_edges,
    const float *host_to_edges,
    GPUNeighborGraph &graph,
    CategoryProfiler *profiler = nullptr);

void gpu_preprocess_free(GPUNeighborGraph &graph);

#endif // GPU_PREPROCESS_HPP
