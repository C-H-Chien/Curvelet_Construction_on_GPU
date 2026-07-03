#ifndef GPU_PREPROCESS_HPP
#define GPU_PREPROCESS_HPP

#include <cstddef>

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

//> Neighbor index storage layout after preprocessing.
enum class GPUNeighborLayout {
    CSR,        //> compact neighbor_offsets + neighbor_ids[]
    FixedRow    //> edgeLookList-style: neighbor_list[anchor * slots + k], unused slots = -1
};

//> How to build the CSR neighbor graph (single-pass is default).
enum class GPUNeighborCSRStrategy {
    SinglePass,     //> discover once -> stage -> scan -> compact
    TwoPass         //> count kernel + fill kernel (rediscover neighbors)
};

//> How to build the fixed-row neighbor table.
enum class GPUFixedRowBuildStrategy {
    Stage,  //> one thread/anchor discovers via discover_and_stage_neighbors_kernel
    Warp    //> one warp/anchor discovers cells in parallel, writes fixed-row directly
};

struct GPUPreprocessConfig {
    int device_id = 0;
    int num_edges = 0;
    int sz_edge_data = 4;
    unsigned neighbor_radius = 3;   // 7x7 window, matches CPU nr = 3
    float rad = 3.5f;
    bool copy_to_host = false;
    GPUNeighborLayout neighbor_layout = GPUNeighborLayout::CSR;
    GPUNeighborCSRStrategy csr_strategy = GPUNeighborCSRStrategy::SinglePass;
    //> Max neighbors staged per anchor (7x7 window, ~2 edges/pixel -> 98)
    unsigned max_candidates = 64;
    //> Two-pass count_neighbors_kernel: threads per block (each thread handles one anchor).
    int neighbor_count_threads = 1;
    //> Two-pass fill_csr_neighbors_kernel: threads per block (each thread handles one anchor).
    int neighbor_fill_threads = 1;
    //> Single-pass stage/compact kernels: threads per block (each thread handles one anchor).
    int neighbor_stage_threads = 1;
    //> Fixed-row warp build: warps per block in discover_fixed_row_warp_kernel (each warp = one anchor).
    int neighbor_warps_per_block = 1;
    GPUFixedRowBuildStrategy fixed_row_build = GPUFixedRowBuildStrategy::Warp;
};

//> Neighbor graph after preprocessing (CSR or fixed-row layout).
struct GPUNeighborGraph {
    GPUNeighborLayout layout = GPUNeighborLayout::CSR;

    int num_edges                = 0;   //> number of edges
    int total_neighbor_pairs     = 0;   //> CSR: neighbor_offsets[num_edges]; sum of neighbor counts
    unsigned max_neighbor_degree = 0;

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

    //> (Optional) host mirrors when copy_to_host = true
    float *host_edges           = nullptr;
    int *host_neighbor_offsets  = nullptr;
    int *host_neighbor_ids      = nullptr;
    float *host_neighbor_dist2  = nullptr;
    int *host_neighbor_list       = nullptr;
    float *host_neighbor_dist2_row = nullptr;
    int *host_neighbor_counts     = nullptr;
};

struct GPUPreprocessResult {
    //> CSR layout
    GPUNeighborGraph csr;
};

//> Build neighbor data on the GPU from a host-side third-order edge array.
//> Optional profiler aggregates timings across pipeline phases (spatial index, neighbor graph, etc.).
bool gpu_preprocess_build(
    const GPUPreprocessConfig &cfg,
    const float *host_to_edges,
    GPUPreprocessResult &result,
    CategoryProfiler *profiler = nullptr);

void gpu_preprocess_free(GPUPreprocessResult &result);

#endif // GPU_PREPROCESS_HPP
