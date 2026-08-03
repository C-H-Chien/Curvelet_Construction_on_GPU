#include "gpu_edge_chain_growth.hpp"
#include "gpu_common.hpp"
#include "timer.hpp"
#include "indices.hpp"
#include "device/gpu_edge_chain_growth_kernels.cuh"

#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

namespace {

void profile_lap(CategoryProfiler *profiler, TimerCategory cat, const char *detail)
{
    if (profiler != nullptr) {
        profiler->lap(cat, detail);
    }
}

size_t warp_smem_floats_per_warp(int smem_mode, int bundle_cells)
{
    if (smem_mode == kWarpSmemLaneWs) {
        return static_cast<size_t>(32) * static_cast<size_t>(2) * static_cast<size_t>(bundle_cells);
    }
    return 0;
}

//> Check the largest available dynamic shared memory per thread-block
bool query_max_dynamic_shared_bytes(size_t &max_dyn_bytes)
{
    int device = 0;
    cudacheck(cudaGetDevice(&device));

    int max_optin = 0;
    cudaError_t err = cudaDeviceGetAttribute( &max_optin, cudaDevAttrMaxSharedMemoryPerBlockOptin, device );
    if (err != cudaSuccess || max_optin <= 0) {
        int max_default = 0;
        cudacheck(cudaDeviceGetAttribute( &max_default, cudaDevAttrMaxSharedMemoryPerBlock, device ));
        max_dyn_bytes = static_cast<size_t>(max_default);
        return max_dyn_bytes > 0;
    }
    max_dyn_bytes = static_cast<size_t>(max_optin);
    return true;
}

//> Optionally cache per-lane growth workspaces in shared memory (pairwise bundles stay in global).
//>
//> chain_smem_mode:
//>   "auto" — prefer lane workspace in shared; may reduce warps_per_block
//>   "none" — mode 0; keep requested warps_per_block
//>   "lane" — mode 1; keep requested warps_per_block (fail if too large)
bool choose_warp_smem_config(
    int bundle_cells,
    int requested_warps_per_block,
    const std::string &chain_smem_mode,
    int &smem_mode,
    int &warps_per_block,
    size_t &smem_bytes)
{
    warps_per_block = requested_warps_per_block;
    smem_mode = kWarpSmemNone;
    smem_bytes = 0;

    if (requested_warps_per_block <= 0 || requested_warps_per_block > 32) {
        fprintf(stderr, "choose_warp_smem_config: invalid warps_per_block=%d (use 1..32)\n",
                requested_warps_per_block);
        return false;
    }

    if (chain_smem_mode == "none") {
        return true;
    }

    size_t max_dyn_bytes = 0;
    if (!query_max_dynamic_shared_bytes(max_dyn_bytes)) {
        if (chain_smem_mode == "auto") {
            return true;
        }
        fprintf(stderr, "choose_warp_smem_config: could not query device shared-memory limit\n");
        return false;
    }

    //> Request max dynamic shared memory
    cudacheck(cudaFuncSetAttribute( grow_edge_chains_warp_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(max_dyn_bytes)));

    auto try_mode = [&](int mode, int warps) -> bool {
        const size_t per_warp = warp_smem_floats_per_warp(mode, bundle_cells);
        const size_t bytes = per_warp * static_cast<size_t>(warps) * sizeof(float);
        if (per_warp == 0 || bytes > max_dyn_bytes) {
            return false;
        }
        smem_mode = mode;
        warps_per_block = warps;
        smem_bytes = bytes;
        return true;
    };

    auto fail_manual = [&](int mode) -> bool {
        const size_t per_warp = warp_smem_floats_per_warp(mode, bundle_cells);
        const size_t need = per_warp * static_cast<size_t>(requested_warps_per_block) * sizeof(float);
        fprintf(stderr,
                "choose_warp_smem_config: --chain-smem-mode=%s needs %zu bytes "
                "(%d warps/block) but device max dynamic shared is %zu bytes; "
                "lower --chain-warps-per-block or use --chain-smem-mode none\n",
                chain_smem_mode.c_str(), need, requested_warps_per_block, max_dyn_bytes);
        return false;
    };

    if (chain_smem_mode == "lane") {
        if (try_mode(kWarpSmemLaneWs, requested_warps_per_block)) {
            return true;
        }
        return fail_manual(kWarpSmemLaneWs);
    }

    //> auto: lane workspace in shared, else global (may reduce warps)
    if (try_mode(kWarpSmemLaneWs, requested_warps_per_block)) {
        return true;
    }
    for (int w = requested_warps_per_block - 1; w >= 1; --w) {
        if (try_mode(kWarpSmemLaneWs, w)) {
            fprintf(stderr,
                    "grow_edge_chains_warp: lane workspace shared cache; reduced warps_per_block %d -> %d\n",
                    requested_warps_per_block, w);
            return true;
        }
    }
    fprintf(stderr,
            "grow_edge_chains_warp: shared memory too small; using global scratch for lane workspace\n");
    smem_mode = kWarpSmemNone;
    warps_per_block = requested_warps_per_block;
    smem_bytes = 0;
    return true;
}

void compute_scratch_sizes(
    int smem_mode,
    int slots_per_anchor,
    int bundle_cells,
    int chain_width,
    size_t &floats_per_anchor,
    size_t &uints_per_anchor)
{
    //> Per-seed final working min/max for forward and backward growth directions
    floats_per_anchor = static_cast<size_t>(4) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(bundle_cells);
    //> Lane workspaces stay in global memory only when shared caching is unavailable
    if (smem_mode == kWarpSmemNone) {
        floats_per_anchor += static_cast<size_t>(32) * static_cast<size_t>(2) * static_cast<size_t>(bundle_cells);
    }
    //> candidate chains for both directions + dedup table
    uints_per_anchor = static_cast<size_t>(3) * static_cast<size_t>(slots_per_anchor) * static_cast<size_t>(chain_width);
}

void launch_grow_edge_chains_warp(
    const CurveletParams &params,
    const GPUNeighborGraph &graph,
    const GPUCurveBundleStorage &bundles,
    GPUCurveletChainStorage &storage)
{
    //> TODO: Remove storage.warp_wraps_per_block; just use params.chain_warps_per_block
    const int warps_per_block = storage.warp_warps_per_block > 0 ? (storage.warp_warps_per_block) : (params.chain_warps_per_block);
    if (warps_per_block <= 0 || warps_per_block > 32) {
        fprintf(stderr, "launch_grow_edge_chains_warp: invalid chain_warps_per_block=%d (use 1..32)\n", warps_per_block);
        return;
    }

    //> Phase 1: warp-per-anchor growth for both forward and backward into dual candidate buffers
    const int grow_threads_per_block = warps_per_block * 32;
    const int grow_num_blocks = div_up(graph.num_edges, warps_per_block);

    grow_edge_chains_warp_kernel<<<grow_num_blocks, grow_threads_per_block, storage.warp_smem_bytes>>>(
        graph.num_edges,
        storage.slots_per_anchor,
        storage.bundle_cells,
        storage.group_max_sz,
        storage.chain_width,
        params.edge_data_sz,
        warps_per_block,
        storage.warp_smem_mode,
        storage.scratch_floats_per_anchor,
        storage.scratch_uints_per_anchor,
        graph.dev_edges,
        graph.dev_neighbor_list,
        graph.dev_neighbor_counts,
        bundles.dev_bundle_min_ks,
        bundles.dev_bundle_max_ks,
        bundles.dev_hyp_look_edge,
        storage.dev_scratch_f,
        storage.dev_scratch_u);
    cudacheck(cudaGetLastError());
}

void launch_dedup_record_edge_chains(
    const CurveletParams &params,
    const GPUNeighborGraph &graph,
    const GPUCurveBundleStorage &bundles,
    GPUCurveletChainStorage &storage)
{
    //> Phase 2: one thread per anchor, deduplicate and record accepted curvelets
    const int dedup_threads_per_block = params.dedup_threads_per_block;
    const int dedup_num_blocks = div_up(graph.num_edges, dedup_threads_per_block);

    dedup_record_edge_chains_kernel<<<dedup_num_blocks, dedup_threads_per_block>>>(
        graph.num_edges,
        storage.slots_per_anchor,
        storage.bundle_cells,
        bundles.curves_num_in_bundle_pixel,
        bundles.curves_num_in_bundle_theta,
        storage.group_max_sz,
        storage.max_per_anchor,
        storage.chain_width,
        static_cast<float>(params.sx),
        static_cast<float>(params.st),
        storage.scratch_floats_per_anchor,
        storage.scratch_uints_per_anchor,
        graph.dev_neighbor_counts,
        storage.dev_scratch_f,
        storage.dev_scratch_u,
        storage.dev_edge_chain_final,
        storage.dev_anchor_chain_count,
        storage.dev_curvelet_info);
    cudacheck(cudaGetLastError());
}

} // namespace

bool gpu_allocate_edge_chains(
    const GPUNeighborGraph &graph,
    const GPUCurveBundleStorage &bundles,
    int group_max_sz,
    int chain_warps_per_block,
    const std::string &chain_smem_mode,
    GPUCurveletChainStorage &storage)
{
    //> Sanity check: make sure that the pairwise bundles exist
    if (bundles.dev_bundle_min_ks == nullptr || bundles.dev_hyp_look_edge == nullptr) {
        fprintf(stderr, "gpu_allocate_edge_chains: pairwise bundles not allocated\n");
        return false;
    }
    if (group_max_sz <= 0) {
        fprintf(stderr, "gpu_allocate_edge_chains: invalid group_max_sz=%d\n", group_max_sz);
        return false;
    }

    storage.num_edges = graph.num_edges;
    storage.slots_per_anchor = graph.neighbor_slots_per_anchor;
    storage.group_max_sz = group_max_sz;
    storage.chain_width = group_max_sz + 1;
    storage.max_per_anchor = (storage.slots_per_anchor + 1) * 2;
    storage.max_curvelets = storage.num_edges * storage.max_per_anchor;
    storage.bundle_cells = bundles.bundle_cells;

    //> Choose the warp-based shared memory configuration
    const int req_warps = (chain_warps_per_block > 0) ? chain_warps_per_block : 4;
    if (!choose_warp_smem_config(storage.bundle_cells, req_warps,
                                 chain_smem_mode,
                                 storage.warp_smem_mode, storage.warp_warps_per_block,
                                 storage.warp_smem_bytes)) {
        return false;
    }
#if VERBOSE
    std::cout << "chain_smem_mode (request): " << chain_smem_mode << std::endl;
    std::cout << "warp_smem_mode:            " << storage.warp_smem_mode << std::endl;
    std::cout << "warp_warps_per_block:      " << storage.warp_warps_per_block << std::endl;
    std::cout << "warp_smem_bytes:           " << storage.warp_smem_bytes << std::endl;
#endif

    size_t floats_per_anchor = 0;
    size_t uints_per_anchor = 0;
    compute_scratch_sizes(storage.warp_smem_mode, storage.slots_per_anchor,
                          storage.bundle_cells, storage.chain_width,
                          floats_per_anchor, uints_per_anchor);

    //> Check if the scratch buffers need to be reallocated
    const bool need_realloc =
        (storage.dev_scratch_f == nullptr) ||
        (storage.scratch_floats_per_anchor != floats_per_anchor) ||
        (storage.scratch_uints_per_anchor != uints_per_anchor);

    if (need_realloc && storage.dev_scratch_f != nullptr) {
        cudaFree(storage.dev_scratch_f);        
        cudaFree(storage.dev_scratch_u);
        storage.dev_scratch_f = nullptr;
        storage.dev_scratch_u = nullptr;
    }

    storage.scratch_floats_per_anchor = floats_per_anchor;
    storage.scratch_uints_per_anchor = uints_per_anchor;

    if (storage.dev_edge_chain_final == nullptr) {
        const size_t chain_words = static_cast<size_t>(storage.max_curvelets) * static_cast<size_t>(storage.chain_width);
        const size_t info_words = static_cast<size_t>(storage.max_curvelets) * static_cast<size_t>(GPU_CURVELET_INFO_WIDTH);

        cudacheck(cudaMalloc(&storage.dev_edge_chain_final, chain_words * sizeof(unsigned)));
        cudacheck(cudaMalloc(&storage.dev_anchor_chain_count, static_cast<size_t>(storage.num_edges) * sizeof(unsigned)));
        cudacheck(cudaMalloc(&storage.dev_curvelet_info, info_words * sizeof(float)));

        cudacheck(cudaMemset(storage.dev_edge_chain_final, 0, chain_words * sizeof(unsigned)));
        cudacheck(cudaMemset(storage.dev_anchor_chain_count, 0, static_cast<size_t>(storage.num_edges) * sizeof(unsigned)));
        cudacheck(cudaMemset(storage.dev_curvelet_info, 0, info_words * sizeof(float)));
    }

    if (storage.dev_scratch_f == nullptr) {
        const size_t scratch_f_words = static_cast<size_t>(storage.num_edges) * storage.scratch_floats_per_anchor;
        const size_t scratch_u_words = static_cast<size_t>(storage.num_edges) * storage.scratch_uints_per_anchor;
        cudacheck(cudaMalloc(&storage.dev_scratch_f, scratch_f_words * sizeof(float)));
        cudacheck(cudaMalloc(&storage.dev_scratch_u, scratch_u_words * sizeof(unsigned)));
    }

    return true;
}

void gpu_curvelet_free_chains(GPUCurveletChainStorage &storage)
{
    if (storage.dev_edge_chain_final != nullptr) {
        cudaFree(storage.dev_edge_chain_final);
        storage.dev_edge_chain_final = nullptr;
    }
    if (storage.dev_anchor_chain_count != nullptr) {
        cudaFree(storage.dev_anchor_chain_count);
        storage.dev_anchor_chain_count = nullptr;
    }
    if (storage.dev_curvelet_info != nullptr) {
        cudaFree(storage.dev_curvelet_info);
        storage.dev_curvelet_info = nullptr;
    }
    if (storage.dev_scratch_f != nullptr) {
        cudaFree(storage.dev_scratch_f);
        storage.dev_scratch_f = nullptr;
    }
    if (storage.dev_scratch_u != nullptr) {
        cudaFree(storage.dev_scratch_u);
        storage.dev_scratch_u = nullptr;
    }
}

bool gpu_grow_edge_chains_main(
    const CurveletParams &params,
    const GPUNeighborGraph &graph,
    const GPUCurveBundleStorage &bundles,
    GPUCurveletChainStorage &storage,
    GPUCurveletChainResult &result,
    CategoryProfiler *profiler)
{
    //> Sanity check: the neighbor graph must follow the fixed-row layout
    if (graph.layout != GPUNeighborLayout::FixedRow) {
        fprintf(stderr, "gpu_grow_edge_chains_main: requires fixed-row neighbor layout\n");
        return false;
    }

    //> Allocate edge chain storage if not already done
    const bool first_alloc = (storage.dev_edge_chain_final == nullptr);
    if ( !gpu_allocate_edge_chains(graph, bundles, static_cast<int>(params.group_max_sz),
                                    params.chain_warps_per_block,
                                    params.chain_smem_mode, storage) ) {
        return false;
    }
    profile_lap(profiler, TimerCategory::MemoryAlloc, first_alloc ? "edge-chain output and scratch buffers" : "edge-chain buffers reused");

    //> Phase 1: warp-per-anchor growth into candidate / working-bundle buffers
    launch_grow_edge_chains_warp(params, graph, bundles, storage);
    cudacheck(cudaDeviceSynchronize());
    profile_lap(profiler, TimerCategory::Kernel, "grow_edge_chains_warp (phase 1)");

    //> Phase 2: one thread per anchor — dedup and record accepted curvelets
    launch_dedup_record_edge_chains(params, graph, bundles, storage);
    cudacheck(cudaDeviceSynchronize());
    profile_lap(profiler, TimerCategory::Kernel, "dedup_record_edge_chains (phase 2)");

    //> Download the anchor chain counts from the device to the host
    std::vector<unsigned> host_counts(static_cast<size_t>(storage.num_edges));
    cudacheck(cudaMemcpy(host_counts.data(), storage.dev_anchor_chain_count, static_cast<size_t>(storage.num_edges) * sizeof(unsigned), cudaMemcpyDeviceToHost));
    profile_lap(profiler, TimerCategory::DataTransfer, "anchor chain counts D->H");

    //> Compute the total number of curvelets
    unsigned total = 0;
    for (unsigned c : host_counts) {
        total += c;
    }
    result.num_curvelets = total;

    return true;
}

bool gpu_download_compact_chains(
    const GPUCurveletChainStorage &storage,
    std::vector<int> &host_chains,
    std::vector<float> &host_info,
    unsigned &num_curvelets)
{
    if (storage.dev_edge_chain_final == nullptr || storage.dev_anchor_chain_count == nullptr) {
        fprintf(stderr, "gpu_download_compact_chains: chain storage not allocated\n");
        return false;
    }

    std::vector<unsigned> host_counts(static_cast<size_t>(storage.num_edges));
    cudacheck(cudaMemcpy(host_counts.data(), storage.dev_anchor_chain_count,
                         static_cast<size_t>(storage.num_edges) * sizeof(unsigned),
                         cudaMemcpyDeviceToHost));

    unsigned total = 0;
    for (unsigned c : host_counts) {
        total += c;
    }
    num_curvelets = total;

    const size_t chain_words = static_cast<size_t>(storage.max_curvelets) * static_cast<size_t>(storage.chain_width);
    const size_t info_words = static_cast<size_t>(storage.max_curvelets) * static_cast<size_t>(GPU_CURVELET_INFO_WIDTH);

    std::vector<unsigned> packed_chains(chain_words);
    std::vector<float> packed_info(info_words);
    cudacheck(cudaMemcpy(packed_chains.data(), storage.dev_edge_chain_final,
                         chain_words * sizeof(unsigned), cudaMemcpyDeviceToHost));
    cudacheck(cudaMemcpy(packed_info.data(), storage.dev_curvelet_info,
                         info_words * sizeof(float), cudaMemcpyDeviceToHost));

    host_chains.assign(static_cast<size_t>(total) * static_cast<size_t>(storage.chain_width), 0);
    host_info.assign(static_cast<size_t>(total) * static_cast<size_t>(GPU_CURVELET_INFO_WIDTH), 0.f);

    unsigned dst = 0;
    for (int te = 0; te < storage.num_edges; te++) {
        const unsigned n = host_counts[static_cast<size_t>(te)];
        for (unsigned k = 0; k < n; k++) {
            const unsigned src = static_cast<unsigned>(te) * static_cast<unsigned>(storage.max_per_anchor) + k;
            for (int j = 0; j < storage.chain_width; j++) {
                host_chains[static_cast<size_t>(dst) * storage.chain_width + static_cast<size_t>(j)] =
                    static_cast<int>(packed_chains[static_cast<size_t>(src) * storage.chain_width + static_cast<size_t>(j)]);
            }
            for (int col = 0; col < GPU_CURVELET_INFO_WIDTH; col++) {
                host_info[static_cast<size_t>(dst) * GPU_CURVELET_INFO_WIDTH + static_cast<size_t>(col)] =
                    packed_info[static_cast<size_t>(src) * GPU_CURVELET_INFO_WIDTH + static_cast<size_t>(col)];
            }
            dst++;
        }
    }
    return true;
}
