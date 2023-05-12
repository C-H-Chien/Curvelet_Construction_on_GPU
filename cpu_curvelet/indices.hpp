#ifndef INDICES_HPP
#define INDICES_HPP
// macros for flexible axis

// Enablers
#define ReadThirdOrderEdgesFromFile         (1)

// Third-Order Edge List for Pre-Processing
#define TO_edges(i, j)                              TO_edges[(i) * sz_edge_data + (j)]
#define edgeLookList(i, j)                      edgeLookList[(i) * ((_sz_edge_data+1)*32) + (j)]
//#define unsorted_edgeLookList(i, j)    unsorted_edgeLookList[(i) * ((_sz_edge_data+1)*32) + (j)]

// Build Curvelets
#define _edgeLookList(i, j)             _edgeLookList[(i) * ((_sz_edge_data+1)*(_max_num_look_edges+1)) + (j)]
#define bundle_min_ks(i, j)             bundle_min_ks[(i) * (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) + (j)]
#define bundle_max_ks(i, j)             bundle_max_ks[(i) * (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) + (j)]
#define cmp_bundle_min_ks(i, j)     cmp_bundle_min_ks[(i) * (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) + (j)]
#define cmp_bundle_max_ks(i, j)     cmp_bundle_max_ks[(i) * (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) + (j)]

// Edge Chains
#define edge_chain_target(i, j)     edge_chain_target[(i) * (_group_max_sz+1) + (j)]
#define edge_chain_final(i, j)       edge_chain_final[(i) * (_group_max_sz+1) + (j)]


#endif // INDICES_HPP