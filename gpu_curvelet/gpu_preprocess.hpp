#ifndef GPU_PREPROCESS_HPP
#define GPU_PREPROCESS_HPP

#include <cstddef>

// GPU preprocessing for the neighbor-search stage of curvelet construction.
//
// Two types of output layouts are explored (select via GPUPreprocessConfig::output_layout):
//
// 1) EdgeLookList (legacy / CPU-compatible padded table)
//    stride = (sz_edge_data + 1) * look_slots
//    row i, slot 0:           [id, x, y, orientation, strength] of anchor edge i
//    row i, slot j (j >= 1):  [id, x, y, orientation, strength] of j-th nearest neighbor
//    unused slots have id = -1
//
// 2) NeighborCSR (compressed sparse rows (CSR) format)
//    Global edge attributes live once in dev_edges (num_edges x sz_edge_data).
//    Per-anchor neighbors are stored in CSR form:
//      neighbor_offsets[anchor] .. neighbor_offsets[anchor+1]  -> slice in neighbor_ids[]
//      neighbor_ids[k] is sorted by ascending squared distance to its anchor.
//      neighbor_dist2[k] holds that squared distance (optional but useful).
//
// Internally both paths share the same spatial cell index (also CSR over pixel cells).

enum class GPUPreprocessOutput {
    EdgeLookList,   //> padded table only
    NeighborCSR,    //> CSR format of compact adjacency only
    Both            //> build both for comparison / validation
};

//> How to build the CSR neighbor graph (single-pass is default).
enum class GPUNeighborCSRStrategy {
    SinglePass,     //> discover once → stage → scan → compact
    TwoPass,        //> count kernel + fill kernel (rediscover neighbors)
    CompareCSR      //> build both CSR strategies and compare
};

struct GPUPreprocessConfig {
    int device_id = 0;
    int num_edges = 0;
    int sz_edge_data = 4;
    unsigned look_slots = 64;
    unsigned neighbor_radius = 3;   // 7x7 window, matches CPU nr = 3
    float rad = 3.5f;
    bool copy_to_host = true;
    GPUPreprocessOutput output_layout = GPUPreprocessOutput::NeighborCSR;
    GPUNeighborCSRStrategy csr_strategy = GPUNeighborCSRStrategy::SinglePass;
    //> Max neighbors staged per anchor (7x7 window, ~2 edges/pixel → 98; 128 is a safe default).
    //> Does not need to be a power of two.
    unsigned max_candidates = 128;
};

//> Compact per-anchor neighbor graph
struct GPUNeighborGraph {
    int num_edges                = 0;   //> number of edges
    int total_neighbor_pairs     = 0;   //> neighbor_offsets[num_edges]
    unsigned max_neighbor_degree = 0;

    //> Global edge table on device (row-major: x, y, orientation, strength)
    float *dev_edges = nullptr;

    //> neighbors of anchor a are stored (based on the distance to that anchor) in the following way:
    //> neighbor_ids[ neighbor_offsets[a] .. neighbor_offsets[a+1] )
    int *dev_neighbor_offsets   = nullptr;    //> size of num_edges + 1
    int *dev_neighbor_ids       = nullptr;    //> size of total_neighbor_pairs
    float *dev_neighbor_dist2   = nullptr;    //> size of total_neighbor_pairs

    //> (Optional) mirroring the host when copy_to_host = true
    float *host_edges           = nullptr;
    int *host_neighbor_offsets  = nullptr;
    int *host_neighbor_ids      = nullptr;
    float *host_neighbor_dist2  = nullptr;
};

struct GPUPreprocessResult {
    //> EdgeLookList layout
    float *dev_edgeLookList = nullptr;
    float *host_edgeLookList = nullptr;
    unsigned max_look_edge_num = 0;
    int edge_look_stride = 0;

    //> CSR layout
    GPUNeighborGraph csr;
};

//> Build neighbor data on the GPU from a host-side third-order edge array.
bool gpu_preprocess_build( const GPUPreprocessConfig &cfg, const float *host_to_edges, GPUPreprocessResult &result);

void gpu_preprocess_free(GPUPreprocessResult &result);

//> Compare CSR neighbors against edgeLookList slot 1..N for every anchor.
bool gpu_preprocess_compare_csr_to_edge_look_list( const GPUPreprocessResult &result, int num_edges_to_check = -1);

bool gpu_preprocess_compare_csr_graphs(
    const GPUNeighborGraph &a,
    const GPUNeighborGraph &b,
    int num_edges_to_check);

// Legacy entry point (EdgeLookList only).
inline bool gpu_preprocess_build_edge_look_list( const GPUPreprocessConfig &cfg, const float *host_to_edges, GPUPreprocessResult &result )
{
    GPUPreprocessConfig local = cfg;
    local.output_layout = GPUPreprocessOutput::EdgeLookList;
    return gpu_preprocess_build(local, host_to_edges, result);
}

#endif // GPU_PREPROCESS_HPP
