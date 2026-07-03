#ifndef INDICES_HPP
#define INDICES_HPP

// Fixed-row bundle storage: dev_bundle_*[anchor * slots * bundle_cells + slot * bundle_cells + cell]

#define bundle_min_ks_at(buf, slots, cells, anchor, slot, cell) \
    (buf)[((size_t)(anchor) * (size_t)(slots) + (size_t)(slot)) * (size_t)(cells) + (size_t)(cell)]

#define bundle_max_ks_at(buf, slots, cells, anchor, slot, cell) \
    (buf)[((size_t)(anchor) * (size_t)(slots) + (size_t)(slot)) * (size_t)(cells) + (size_t)(cell)]

#define hyp_look_edge_at(buf, slots, anchor, slot) \
    (buf)[(size_t)(anchor) * (size_t)(slots) + (size_t)(slot)]

//> TODO: update these when edge chain growth is implemented
#define cmp_bundle_min_ks(i, j) cmp_bundle_min_ks[(i) * bundle_cells + (j)]
#define cmp_bundle_max_ks(i, j) cmp_bundle_max_ks[(i) * bundle_cells + (j)]

#endif // INDICES_HPP
