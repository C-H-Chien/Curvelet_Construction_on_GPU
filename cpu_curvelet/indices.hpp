#ifndef INDICES_HPP
#define INDICES_HPP
// macros for flexible data access

// Build Curvelets
#define bundle_min_ks(i, j)             bundle_min_ks[(i) * (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) + (j)]
#define bundle_max_ks(i, j)             bundle_max_ks[(i) * (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) + (j)]
#define cmp_bundle_min_ks(i, j)     cmp_bundle_min_ks[(i) * (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) + (j)]
#define cmp_bundle_max_ks(i, j)     cmp_bundle_max_ks[(i) * (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) + (j)]

// Edge Chains
#define edge_chain_target(i, j)     edge_chain_target[(i) * (_group_max_sz+1) + (j)]
#define edge_chain_final(i, j)       _edge_chain_final[(i) * (_group_max_sz+1) + (j)]


#endif // INDICES_HPP
