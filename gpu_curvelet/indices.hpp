#ifndef INDICES_HPP
#define INDICES_HPP

#define VERBOSE (true)

#define bundle_min_ks_at(buf, slots, cells, anchor, slot, cell) \
    (buf)[((size_t)(anchor) * (size_t)(slots) + (size_t)(slot)) * (size_t)(cells) + (size_t)(cell)]

#define bundle_max_ks_at(buf, slots, cells, anchor, slot, cell) \
    (buf)[((size_t)(anchor) * (size_t)(slots) + (size_t)(slot)) * (size_t)(cells) + (size_t)(cell)]

#define hyp_look_edge_at(buf, slots, anchor, slot) \
    (buf)[(size_t)(anchor) * (size_t)(slots) + (size_t)(slot)]

//> Compare buffers used during chain growth (device scratch; working min/max only).
#define cmp_bundle_min_ks(i, j) cmp_bundle_min_ks[(i) * bundle_cells + (j)]
#define cmp_bundle_max_ks(i, j) cmp_bundle_max_ks[(i) * bundle_cells + (j)]

#define edge_chain_target(i, j) edge_chain_target[(i) * (group_max_sz + 1) + (j)]
#define edge_chain_final(row, j) dev_edge_chain_final[(size_t)(row) * (size_t)(chain_width) + (size_t)(j)]

#endif // INDICES_HPP
