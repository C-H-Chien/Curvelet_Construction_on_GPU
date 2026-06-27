#ifndef PREPROCESS_CSR_HPP
#define PREPROCESS_CSR_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "curvelet_utils.hpp"

// CPU neighbor graph in CSR form (mirrors GPUNeighborGraph on the host).
template<typename T>
struct CPUNeighborGraph {
    int num_edges = 0;
    int total_neighbor_pairs = 0;
    unsigned max_num_of_neighbors = 0;

    // Global edge table (row-major: x, y, orientation, strength)
    std::vector<T> edges;

    // neighbor_ids[ neighbor_offsets[a] .. neighbor_offsets[a+1] ), sorted by distance
    std::vector<int> neighbor_offsets;
    std::vector<int> neighbor_ids;
    std::vector<T> neighbor_dist2;
};

struct SpatialIndexCSR {
    std::vector<std::int64_t> unique_cells;
    std::vector<int> cell_starts;
    std::vector<int> cell_counts;
    std::vector<int> sorted_edge_ids;
    int num_cells = 0;
};

inline std::int64_t pack_cell_key(int x, int y)
{
    return (static_cast<std::int64_t>(y) << 32) |
           (static_cast<std::uint32_t>(x));
}

inline int lower_bound_cell_keys(const std::int64_t *keys, int n, std::int64_t key)
{
    int lo = 0;
    int hi = n;
    while (lo < hi) {
        const int mid = lo + ((hi - lo) >> 1);
        if (keys[mid] < key) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

template<typename T>
void build_spatial_index_cpu(
    int num_edges,
    int sz_edge_data,
    const T *edges,
    SpatialIndexCSR &index)
{
    index.unique_cells.clear();
    index.cell_starts.clear();
    index.cell_counts.clear();
    index.sorted_edge_ids.clear();
    index.num_cells = 0;

    if (num_edges <= 0) {
        return;
    }

    std::vector<std::int64_t> cell_keys(static_cast<size_t>(num_edges));
    std::vector<int> perm(static_cast<size_t>(num_edges));
    for (int i = 0; i < num_edges; ++i) {
        perm[static_cast<size_t>(i)] = i;
        const T x = edges[i * sz_edge_data + 0];
        const T y = edges[i * sz_edge_data + 1];
        const int cx = static_cast<int>(std::lrint(x));
        const int cy = static_cast<int>(std::lrint(y));
        cell_keys[static_cast<size_t>(i)] = pack_cell_key(cx, cy);
    }

    std::sort(perm.begin(), perm.end(), [&](int a, int b) {
        return cell_keys[static_cast<size_t>(a)] < cell_keys[static_cast<size_t>(b)];
    });

    index.sorted_edge_ids.resize(static_cast<size_t>(num_edges));
    for (int p = 0; p < num_edges; ++p) {
        index.sorted_edge_ids[static_cast<size_t>(p)] = perm[static_cast<size_t>(p)];
    }

    int p = 0;
    while (p < num_edges) {
        const std::int64_t key = cell_keys[static_cast<size_t>(perm[static_cast<size_t>(p)])];
        const int start = p;
        while (p < num_edges &&
               cell_keys[static_cast<size_t>(perm[static_cast<size_t>(p)])] == key) {
            ++p;
        }
        index.unique_cells.push_back(key);
        index.cell_counts.push_back(p - start);
        index.cell_starts.push_back(start);
    }
    index.num_cells = static_cast<int>(index.unique_cells.size());
}

template<typename T>
void discover_neighbors_csr(
    int te_idx,
    const T *edges,
    int sz_edge_data,
    const SpatialIndexCSR &spatial,
    T rad_sqr,
    unsigned neighbor_radius,
    unsigned max_candidates,
    std::vector<int> &out_ids,
    std::vector<T> &out_dist2)
{
    out_ids.clear();
    out_dist2.clear();

    const T te_x = edges[te_idx * sz_edge_data + 0];
    const T te_y = edges[te_idx * sz_edge_data + 1];
    const int cx = static_cast<int>(std::lrint(te_x));
    const int cy = static_cast<int>(std::lrint(te_y));

    struct Candidate {
        T dist;
        int id;
    };
    std::vector<Candidate> candidates;
    candidates.reserve(max_candidates);

    for (int py = cy - static_cast<int>(neighbor_radius); py <= cy + static_cast<int>(neighbor_radius); ++py) {
        for (int px = cx - static_cast<int>(neighbor_radius); px <= cx + static_cast<int>(neighbor_radius); ++px) {

            const std::int64_t key = pack_cell_key( px, py );
            const int cell_idx = lower_bound_cell_keys( spatial.unique_cells.data(), spatial.num_cells, key );
            if (cell_idx >= spatial.num_cells || spatial.unique_cells[static_cast<size_t>(cell_idx)] != key) {
                continue;
            }

            const int start = spatial.cell_starts[static_cast<size_t>(cell_idx)];
            const int end = start + spatial.cell_counts[static_cast<size_t>(cell_idx)];
            for (int k = start; k < end; ++k) {
                const int ne_idx = spatial.sorted_edge_ids[static_cast<size_t>(k)];
                if (ne_idx == te_idx) {
                    continue;
                }

                const T ne_x = edges[ne_idx * sz_edge_data + 0];
                const T ne_y = edges[ne_idx * sz_edge_data + 1];
                const T dist = sq_dist<T>(te_x, te_y, ne_x, ne_y);
                if (dist > rad_sqr) {
                    continue;
                }

                if (static_cast<unsigned>(candidates.size()) < max_candidates) {
                    candidates.push_back({dist, ne_idx});
                }
            }
        }
    }

    std::sort(candidates.begin(), candidates.end(),
        [](const Candidate &a, const Candidate &b) {
            return a.dist < b.dist;
        }
    );

    out_ids.reserve(candidates.size());
    out_dist2.reserve(candidates.size());
    for (const Candidate &c : candidates) {
        out_ids.push_back(c.id);
        out_dist2.push_back(c.dist);
    }
}

template<typename T>
bool build_neighbor_csr_graph(
    int num_edges,
    int sz_edge_data,
    const T *edges,
    T rad,
    unsigned neighbor_radius,
    unsigned max_candidates,
    int nthreads,
    CPUNeighborGraph<T> &graph)
{
    if (num_edges <= 0 || sz_edge_data < 4 || edges == nullptr) {
        return false;
    }

#ifdef _OPENMP
    if (nthreads > 0) {
        omp_set_num_threads(nthreads);
    }
#endif

    //> Map the third-order edge locations to cells in a spatial index
    SpatialIndexCSR spatial;
    build_spatial_index_cpu(num_edges, sz_edge_data, edges, spatial);

    std::vector<std::vector<int>> staged_ids(static_cast<size_t>(num_edges));
    std::vector<std::vector<T>> staged_dist2(static_cast<size_t>(num_edges));
    unsigned truncated_anchors = 0;

    //> Identify the neighbors for each third-order edge 
    //> This should give (i) a data array and (ii) an offset array
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(nthreads > 1)
#endif
    for (int te = 0; te < num_edges; ++te) {
        discover_neighbors_csr( te, edges, sz_edge_data, spatial, rad * rad, neighbor_radius, max_candidates,
                                staged_ids[static_cast<size_t>(te)], staged_dist2[static_cast<size_t>(te)]);

        if (staged_ids[static_cast<size_t>(te)].size() >= max_candidates) {
#ifdef _OPENMP
            #pragma omp atomic
#endif
            ++truncated_anchors;
        }
    }

    //> Initialize the CSR structure for neighbor graph
    graph.num_edges = num_edges;
    graph.edges.assign(edges, edges + static_cast<size_t>(num_edges) * sz_edge_data);
    graph.neighbor_offsets.assign(static_cast<size_t>(num_edges) + 1, 0);
    graph.max_num_of_neighbors = 0;

    //> Count the number of neighbors for each third-order edge
    //> This completes the offset array
    for (int te = 0; te < num_edges; ++te) {
        const unsigned deg = static_cast<unsigned>(staged_ids[static_cast<size_t>(te)].size());
        graph.neighbor_offsets[static_cast<size_t>(te) + 1] = graph.neighbor_offsets[static_cast<size_t>(te)] + static_cast<int>(staged_ids[static_cast<size_t>(te)].size());
        if (deg > graph.max_num_of_neighbors) {
            graph.max_num_of_neighbors = deg;
        }
    }

    graph.total_neighbor_pairs = graph.neighbor_offsets[static_cast<size_t>(num_edges)];
    graph.neighbor_ids.resize(static_cast<size_t>(graph.total_neighbor_pairs));
    graph.neighbor_dist2.resize(static_cast<size_t>(graph.total_neighbor_pairs));

    //> Fill in the CSR structure for neighbor graph data arrays
    //> This completes the data array
    for (int te = 0; te < num_edges; ++te) {
        const int base = graph.neighbor_offsets[static_cast<size_t>(te)];
        const std::vector<int> &ids = staged_ids[static_cast<size_t>(te)];
        const std::vector<T> &d2 = staged_dist2[static_cast<size_t>(te)];
        for (size_t i = 0; i < ids.size(); ++i) {
            graph.neighbor_ids[static_cast<size_t>(base) + i] = ids[i];
            graph.neighbor_dist2[static_cast<size_t>(base) + i] = d2[i];
        }
    }

    if (truncated_anchors > 0) {
        std::cerr << "Warning: " << truncated_anchors
                  << " anchors hit max_candidates=" << max_candidates
                  << "; neighbor lists may be truncated\n";
    }

    std::cout << "CSR neighbor graph: " << num_edges << " anchors, "
              << graph.total_neighbor_pairs << " total pairs, max number of neighbors per anchor: "
              << graph.max_num_of_neighbors << std::endl;
    return true;
}

#endif // PREPROCESS_CSR_HPP
