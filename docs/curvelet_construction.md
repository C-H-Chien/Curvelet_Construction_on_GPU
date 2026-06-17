# Curvelet (Curvel) Construction Workflow

A **curvelet**, or **curvel**, is a short chain of third-order edges that are geometrically consistent: they can be explained by a single smooth, local curve segment with bounded curvature uncertainty. The input is a list of third-order edges (subpixel positions and orientation); the output is a list of edge chains (curvelets).

## Repository layout

| Folder | Description |
|--------|------|
| `original_code/` | Reference C++ from the official released code of the paper |
| `cpu_curvelet/` | Rewritten CPU version (flat arrays) with OpenMP support |
| `gpu_curvelet/` | (NOT COMPLETE) GPU port of the rewrite |
| `test_files/` | Sample third-order edge files |

The rewrite (`cpu_curvelet`) keeps the same core math and replaces object-heavy curve-bundle and curvelet storage (`edgemap`, `curveletmap`, `CC_curve_model_3d`) with flat arrays in the construction kernel — the main step needed for the GPU port in `gpu_curvelet/`. Preprocessing still uses CPU-side maps and heap-allocated `edgel` pointers (`preprocess.hpp`).

## High-level pipeline

```
Third-order edges; for each edge (anchor edge)
        │
        ▼
  Identify its local neighborhood edges
        │
        ▼
  Form pairwise curve bundles between the anchor edge and its neighbor edges
        │
        ▼
  Grow edge chains by cumulative curve bundle intersection
        │
        ▼
  Accept / reject curvelet chains for the anchor edge
```

### Identify Local Neighborhood Edges

Each edge is first hashed into pixel-size cells, enabling efficient lookup of candidates in a 7×7 window. Then, for each edge:

- Compute squared Euclidean distance to each candidate edge in those cells and keep those within radius `nrad`.
- Sort surviving neighbors by ascending distance.

| Role | Function | File |
|------|----------|------|
| Orchestration | `run_curvelet()` | `cpu_curvelet/main.cpp` |
| Spatial hash | `edgeMap` (constructor) | `cpu_curvelet/preprocess.hpp` |
| Neighbor search | `edgeNeighborList::create_edgeLookList()` | `cpu_curvelet/preprocess.hpp` |
| Squared distance | `sq_dist()` | `cpu_curvelet/curvelet_utils.hpp` |

### Form Pairwise Curve Bundles Between the Anchor Edge and its Neighbor Edges

For each edge (anchor edge),
- **Apply a direction filter**: the orientation difference between the anchor edge and its neighboring edges must be below a threshold.
- **Pairwise Circular Arc Hypothesis Formation**: construct curve bundles between pairs of anchor edge and its neighbor edges.
    - For anchor edge `T` and neighbor `L`, all smooth curves passing near both edges (within `dx`, `dt` uncertainty) must have curvature in some interval `[k_min, k_max]` representing a curve bundle.
    - Each neighbor owns a curve bundle with the anchor edge.
    - Geometrically, it is basically saying, given uncertainty in location and orientation of the anchor edge and the neighbor edge, which curvatures at the anchor edge could pass near this neighbor.
    - The determination of the curvature feasibility interval at the anchor edge uses the Circular Arc Bundle Transport (Prosposition 4 of [the paper](https://ieeexplore.ieee.org/abstract/document/8382271)).
- Code:

| Role | Function | File |
|------|----------|------|
| Main loop (driver) | `CurveletCPU::build_curvelets_greedy()` | `cpu_curvelet/cpu_curvelet.hpp` |
| Direction filter | `angle_from_pt_to_pt()`, `dot()` | `cpu_curvelet/curvelet_utils.hpp` |
| Bundle transport | `CurveletCPU::compute_curve_bundle()` | `cpu_curvelet/cpu_curvelet.hpp` |
| Bundle validity | `CurveletCPU::bundle_valid_check()` | `cpu_curvelet/cpu_curvelet.hpp` |
| Read neighbor from table | `retrieve_edge_data_from_edgeLookList()` | `cpu_curvelet/cpu_curvelet.hpp` |

### Grow Edge Chains by Cumulative Curve Bundle Intersection

Start from an anchor <-> neighbor curve bundle, intersect with other anchor <-> neighbor curve bundle incrementally until all curve bundles are chained.
- This operation starts from each anchor <-> neighbor edge pair independently, enumerating many candidate curvelets greedily. Thus, there are multiple intersected, candidate curve bundles per anchor edge.
- Geometrically, the intersection means: curvatures at the anchor edge that are simultaneously consistent with passing near all accepted neighbors.
- Code: 

| Role | Function | File |
|------|----------|------|
| for-loop (driver) | `CurveletCPU::build_curvelets_greedy()` | `cpu_curvelet/cpu_curvelet.hpp` |
| Copy bundle into compare buffer | `move_to_cmp_bundle()` | `cpu_curvelet/cpu_curvelet.hpp` |
| Intersect two bundles | `bundle_intersection()` | `cpu_curvelet/cpu_curvelet.hpp` |
| Test non-empty intersection | `bundle_intersection_valid_check()` | `cpu_curvelet/cpu_curvelet.hpp` |
| Prepend edge (trailing mode) | `chain_push_front()` | `cpu_curvelet/cpu_curvelet.hpp` |

### Finalize Curve Bundles for the Anchor Edge
Accept or reject candidate curve bundles by deduplicating exact repeats as long as the length of the curvelet (number of edges participating the curvelet) is at least 3. 
- Code: 

| Role | Function | File |
|------|----------|------|
| Size + dedup gate | inline in `build_curvelets_greedy()` | `cpu_curvelet/cpu_curvelet.hpp` |
| Duplicate check | `check_curvelet_exist()` | `cpu_curvelet/cpu_curvelet.hpp` |
| Record accepted chain | `record_curvelet_chain()` | `cpu_curvelet/cpu_curvelet.hpp` |
| Compact global output | `compact_curvelet_output()` | `cpu_curvelet/cpu_curvelet.hpp` |
| Write output file | `write_int_array_to_file()` | `cpu_curvelet/main.cpp` |

## Format of the Output Files

`chain.txt`: each output row has `group_max_sz + 1` columns (tab-separated). With default `group_max_sz = 7`, that is 8 columns:

```
[anchor_id+1]  [chain[0]+1]  [chain[1]+1]  ...  [padding 0s]
```
Each number is the ID (1-based) of the third-order edges. For example,
```
1    1    11    10    9    8    7    0
^    ^
|    anchor ID appears again as first chain element
anchor ID 
```

## Input parameters (To be updated)

| Parameter | Meaning |
|-----------|---------|
| `nrad` | Grouping radius around each edge (default 3.5 px) |
| `gap` | Distance between consecutive edges in a link (used in downstream stages only; not in used in `cpu_curvelet`) |
| `dx`, `dt` | Position / orientation uncertainty at the anchor edge |
| `sx`, `st` | Sampling step for the bundle grid (0.1 px, 0.08 rad) |
| `max_k` | Maximum allowed curvature |
| `group_max_sz` | Maximum edges per curvelet (7 in paper; `4` in current `cpu_curvelet/main.cpp`) |
| `curvelet_style` | Direction-pass mode (`2` = anchor-leading bidirectional in current `main.cpp`; `3` = ENO leading + trailing in original) |

## Notes

In the original code, some downstream stages continue past curvelet formation: _(i)_ connect curvelets into longer edge links (`construct_the_link_graph()`), _(ii)_ extract 1D contour fragments (`extract_regular_contours_from_the_link_graph()`), and _(iii)_ fit C¹ polyarcs to chains (`fit_polyarcs_to_all_edgel_chains()`). These are beyond the curvelet construction stage and are not the focus of this repository. See `original_code/form_curvelet_process.cpp` for those functions.
