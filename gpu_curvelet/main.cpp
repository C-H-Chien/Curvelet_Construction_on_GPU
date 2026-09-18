#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <sstream>

#include "data_io.hpp"
#include "param_settings.hpp"
#include "gpu_preprocess.hpp"
#include "gpu_curve_bundle_formation.hpp"
#include "gpu_edge_chain_growth.hpp"
#include "gpu_common.hpp"
#include "timer.hpp"
#include "indices.hpp"

namespace {

struct RunTimingSummary {
    std::string run_id;
    std::string edge_file;
    std::string gpu_name;
    int gpu_id = 0;
    int num_edges = 0;

    std::string neighbor_layout;
    std::string fixed_row_build;
    std::string csr_strategy;
    std::string csr_discover_mode;
    int neighbor_warps_per_block = 0;
    int bundle_warps_per_block = 0;
    int chain_warps_per_block_req = 0;
    int chain_warps_per_block_eff = 0;
    std::string chain_smem_mode_req;
    int chain_smem_mode_eff = -1;
    int dedup_threads_per_block = 0;
    unsigned max_candidates = 0;

    unsigned max_num_of_neighbors = 0;
    int neighbor_slots_per_anchor = 0;
    int total_neighbor_pairs = 0;
    double mean_neighbors = 0.0;
    unsigned valid_pairs = 0;
    unsigned num_curvelets = 0;
    bool ran_bundle_chain = false;

    double preprocess_s = 0.0;
    double preprocess_kernel_s = 0.0;
    double preprocess_mem_s = 0.0;
    double preprocess_xfer_s = 0.0;
    double preprocess_thrust_s = 0.0;

    double bundle_s = 0.0;
    double bundle_kernel_s = 0.0;
    double bundle_mem_s = 0.0;
    double bundle_xfer_s = 0.0;

    double chain_s = 0.0;
    double chain_kernel_s = 0.0;
    double chain_grow_kernel_s = 0.0;   //> grow_edge_chains_warp (phase 1)
    double chain_dedup_kernel_s = 0.0;  //> dedup_record_edge_chains (phase 2)
    double chain_mem_s = 0.0;
    double chain_xfer_s = 0.0;

    double wall_s = 0.0;
};

const char *timing_summary_header()
{
    return
        "run_id,edge_file,gpu_name,gpu_id,num_edges,"
        "neighbor_layout,fixed_row_build,csr_strategy,csr_discover_mode,"
        "neighbor_warps_per_block,bundle_warps_per_block,"
        "chain_warps_per_block_req,chain_warps_per_block_eff,"
        "chain_smem_mode_req,chain_smem_mode_eff,dedup_threads_per_block,max_candidates,"
        "max_num_of_neighbors,neighbor_slots_per_anchor,total_neighbor_pairs,mean_neighbors,"
        "valid_pairs,num_curvelets,ran_bundle_chain,"
        "preprocess_s,preprocess_kernel_s,preprocess_mem_s,preprocess_xfer_s,preprocess_thrust_s,"
        "bundle_s,bundle_kernel_s,bundle_mem_s,bundle_xfer_s,"
        "chain_s,chain_kernel_s,chain_grow_kernel_s,chain_dedup_kernel_s,chain_mem_s,chain_xfer_s,"
        "wall_s";
}

std::string timing_summary_row(const RunTimingSummary &r)
{
    std::ostringstream oss;
    oss << csv_escape(r.run_id) << ','
        << csv_escape(r.edge_file) << ','
        << csv_escape(r.gpu_name) << ','
        << r.gpu_id << ','
        << r.num_edges << ','
        << csv_escape(r.neighbor_layout) << ','
        << csv_escape(r.fixed_row_build) << ','
        << csv_escape(r.csr_strategy) << ','
        << csv_escape(r.csr_discover_mode) << ','
        << r.neighbor_warps_per_block << ','
        << r.bundle_warps_per_block << ','
        << r.chain_warps_per_block_req << ','
        << r.chain_warps_per_block_eff << ','
        << csv_escape(r.chain_smem_mode_req) << ','
        << r.chain_smem_mode_eff << ','
        << r.dedup_threads_per_block << ','
        << r.max_candidates << ','
        << r.max_num_of_neighbors << ','
        << r.neighbor_slots_per_anchor << ','
        << r.total_neighbor_pairs << ','
        << csv_format_seconds(r.mean_neighbors) << ','
        << r.valid_pairs << ','
        << r.num_curvelets << ','
        << (r.ran_bundle_chain ? 1 : 0) << ','
        << csv_format_seconds(r.preprocess_s) << ','
        << csv_format_seconds(r.preprocess_kernel_s) << ','
        << csv_format_seconds(r.preprocess_mem_s) << ','
        << csv_format_seconds(r.preprocess_xfer_s) << ','
        << csv_format_seconds(r.preprocess_thrust_s) << ','
        << csv_format_seconds(r.bundle_s) << ','
        << csv_format_seconds(r.bundle_kernel_s) << ','
        << csv_format_seconds(r.bundle_mem_s) << ','
        << csv_format_seconds(r.bundle_xfer_s) << ','
        << csv_format_seconds(r.chain_s) << ','
        << csv_format_seconds(r.chain_kernel_s) << ','
        << csv_format_seconds(r.chain_grow_kernel_s) << ','
        << csv_format_seconds(r.chain_dedup_kernel_s) << ','
        << csv_format_seconds(r.chain_mem_s) << ','
        << csv_format_seconds(r.chain_xfer_s) << ','
        << csv_format_seconds(r.wall_s);
    return oss.str();
}

bool append_timing_summary_csv(const std::string &path, const RunTimingSummary &row)
{
    if (path.empty()) {
        return true;
    }
    const bool write_header = !file_is_nonempty(path);
    std::ofstream out(path, std::ios::app);
    if (!out) {
        std::cerr << "Error: could not open timing CSV: " << path << std::endl;
        return false;
    }
    if (write_header) {
        out << timing_summary_header() << '\n';
    }
    out << timing_summary_row(row) << '\n';
    return static_cast<bool>(out);
}

std::string make_run_id(const std::string &edge_file)
{
    using Clock = std::chrono::system_clock;
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        Clock::now().time_since_epoch())
                        .count();
    //> Keep basename only for readability in joined tables.
    std::string base = edge_file;
    const auto slash = base.find_last_of("/\\");
    if (slash != std::string::npos) {
        base = base.substr(slash + 1);
    }
    return base + "_" + std::to_string(ms);
}

void fill_stage_totals(const CategoryProfiler &p,
                       double &total_s,
                       double &kernel_s,
                       double &mem_s,
                       double &xfer_s,
                       double *thrust_s = nullptr)
{
    total_s = p.elapsed();
    kernel_s = p.category_total(TimerCategory::Kernel);
    mem_s = p.category_total(TimerCategory::MemoryAlloc);
    xfer_s = p.category_total(TimerCategory::DataTransfer);
    if (thrust_s != nullptr) {
        *thrust_s = p.category_total(TimerCategory::Thrust);
    }
}

double detail_seconds_containing(const CategoryProfiler &p, const char *needle)
{
    if (needle == nullptr || needle[0] == '\0') {
        return 0.0;
    }
    double sum = 0.0;
    for (const auto &r : p.details()) {
        if (r.detail.find(needle) != std::string::npos) {
            sum += r.seconds;
        }
    }
    return sum;
}

} // namespace

bool run_curvelet_gpu(const std::string &out_chain_file, int gpu_id, CurveletParams &params,
                      const std::string &timing_csv, const std::string &timing_detail_csv,
                      const std::string &neighbor_degree_csv, bool neighbor_degree_only)
{
    const unsigned curvelet_style = params.curvelet_style;
    const unsigned out_type = params.out_type;

    const std::string &edge_file = params.edge_file;
    const int edge_data_sz = params.edge_data_sz;

    if (params.chain_smem_mode != "auto" && params.chain_smem_mode != "none" &&
        params.chain_smem_mode != "lane" && params.chain_smem_mode != "bundles" &&
        params.chain_smem_mode != "filter-bundles" && params.chain_smem_mode != "tile-nocut" &&
        params.chain_smem_mode != "tile") {
        std::cerr << "Warning: unknown --chain-smem-mode '" << params.chain_smem_mode
                  << "', using auto (expected: auto/none/lane/bundles/filter-bundles/tile-nocut/tile)\n";
        params.chain_smem_mode = "auto";
    }
    if (params.chain_smem_mode == "tile" || params.chain_smem_mode == "tile-nocut") {
        if (params.chain_tile_workspaces < 1 || params.chain_tile_workspaces > 32) {
            std::cerr << "Error: --chain-tile-workspaces must be 1..32 (got "
                      << params.chain_tile_workspaces << ")\n";
            return false;
        }
    }

    std::cout << "Using scalar type: float (GPU)" << std::endl;

    std::vector<float> TOED_edges;
    if (!read_TO_edges_from_file(edge_file, edge_data_sz, TOED_edges)) {
        return false;
    }
    int edge_num = static_cast<int>(TOED_edges.size() / edge_data_sz);

    cudaDeviceProp prop;
    cudacheck(cudaSetDevice(gpu_id));
    cudaGetDeviceProperties(&prop, gpu_id);
    printf("Device name: %s (Compute capability: %d.%d)\n", prop.name, prop.major, prop.minor);

    RunTimingSummary timing;
    timing.run_id = make_run_id(edge_file);
    timing.edge_file = edge_file;
    timing.gpu_name = prop.name;
    timing.gpu_id = gpu_id;
    timing.num_edges = edge_num;
    timing.neighbor_layout = params.neighbor_layout;
    timing.fixed_row_build = params.fixed_row_build;
    timing.csr_strategy = params.csr_strategy;
    timing.csr_discover_mode = params.csr_discover_mode;
    timing.neighbor_warps_per_block = params.neighbor_warps_per_block;
    timing.bundle_warps_per_block = params.bundle_warps_per_block;
    timing.chain_warps_per_block_req = params.chain_warps_per_block;
    timing.chain_smem_mode_req = params.chain_smem_mode;
    timing.dedup_threads_per_block = params.dedup_threads_per_block;
    timing.max_candidates = params.max_candidates;

    const auto wall_t0 = std::chrono::steady_clock::now();

    //> ================== Preprocess ==================
    CategoryProfiler profiler;
    profiler.set_title("preprocess");
    profiler.start();

    GPUNeighborGraph neighbor_graph;
    if (!gpu_preprocess_build(params, gpu_id, edge_num, TOED_edges.data(), neighbor_graph, &profiler)) {
        return false;
    }
    fill_stage_totals(profiler,
                      timing.preprocess_s,
                      timing.preprocess_kernel_s,
                      timing.preprocess_mem_s,
                      timing.preprocess_xfer_s,
                      &timing.preprocess_thrust_s);
    profiler.summary();

    timing.max_num_of_neighbors = neighbor_graph.max_num_of_neighbors;
    timing.neighbor_slots_per_anchor = neighbor_graph.neighbor_slots_per_anchor;
    timing.total_neighbor_pairs = neighbor_graph.total_neighbor_pairs;
    timing.mean_neighbors = (edge_num > 0) ? static_cast<double>(neighbor_graph.total_neighbor_pairs) / static_cast<double>(edge_num) : 0.0;

    if (!neighbor_degree_csv.empty()) {
        if (neighbor_graph.dev_neighbor_counts == nullptr) {
            std::cerr << "Cannot write --neighbor-degree-csv: neighbor counts not available\n";
            gpu_preprocess_free(neighbor_graph);
            return false;
        }
        std::vector<int> host_counts(static_cast<size_t>(edge_num));
        cudacheck(cudaMemcpy(host_counts.data(), neighbor_graph.dev_neighbor_counts,
                             static_cast<size_t>(edge_num) * sizeof(int), cudaMemcpyDeviceToHost));
        if (!write_neighbor_degree_csv(neighbor_degree_csv, host_counts.data(), edge_num)) {
            gpu_preprocess_free(neighbor_graph);
            return false;
        }
    }

    if (neighbor_degree_only && neighbor_degree_csv.empty()) {
        std::cerr << "Error: --neighbor-degree-only requires --neighbor-degree-csv <file>\n";
        return false;
    }

    if (!timing_detail_csv.empty()) {
        profiler.append_detail_csv(timing_detail_csv, timing.run_id);
    }

    if (neighbor_degree_only) {
        gpu_preprocess_free(neighbor_graph);
        timing.wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_t0).count();
        const std::string summary_line = timing_summary_row(timing);
        std::cout << "\nTIMING_CSV " << summary_line << std::endl;
        if (!append_timing_summary_csv(timing_csv, timing)) {
            return false;
        }
        if (!timing_csv.empty()) {
            std::cout << "Wrote timing summary row to " << timing_csv << std::endl;
        }
        std::cout << "Neighbor-degree-only run complete (bundle/chain skipped)." << std::endl;
        return true;
    }

    if (neighbor_graph.layout == "fixed-row") {
        std::cout << "Fixed-row layout ready: " << neighbor_graph.total_neighbor_pairs
                  << " anchor-neighbor pairs, " << neighbor_graph.neighbor_slots_per_anchor
                  << " slots/anchor, max number of neighbors = " << neighbor_graph.max_num_of_neighbors
                  << std::endl;

        //> ================== Pairwise Curve Bundle Formation ==================
        CategoryProfiler bundle_profiler;
        bundle_profiler.set_title("pairwise_curve_bundles");
        bundle_profiler.start();

        GPUCurveBundleStorage bundle_storage;
        unsigned valid_pairs = 0;
        if (!gpu_form_pairwise_bundles_main(params, neighbor_graph, bundle_storage, valid_pairs, &bundle_profiler)) {
            gpu_curvelet_free_bundles(bundle_storage);
            gpu_preprocess_free(neighbor_graph);
            return false;
        }
        fill_stage_totals(bundle_profiler,
                          timing.bundle_s,
                          timing.bundle_kernel_s,
                          timing.bundle_mem_s,
                          timing.bundle_xfer_s);
        timing.valid_pairs = valid_pairs;
        bundle_profiler.summary();
#if VERBOSE
        std::cout << "Pairwise curve bundles formed (fixed-row warp): " << valid_pairs << " valid pairs" << std::endl;
#endif

        if (!timing_detail_csv.empty()) {
            bundle_profiler.append_detail_csv(timing_detail_csv, timing.run_id);
        }

        //> ================== Edge Chain Growth ==================
        CategoryProfiler chain_profiler;
        chain_profiler.set_title("grow_edge_chains");
        chain_profiler.start();

        GPUCurveletChainStorage chain_storage;
        unsigned num_curvelets = 0;
        if (!gpu_grow_edge_chains_main(params, neighbor_graph, bundle_storage, chain_storage, num_curvelets, &chain_profiler)) {
            gpu_curvelet_free_chains(chain_storage);
            gpu_curvelet_free_bundles(bundle_storage);
            gpu_preprocess_free(neighbor_graph);
            return false;
        }
        fill_stage_totals(chain_profiler,
                          timing.chain_s,
                          timing.chain_kernel_s,
                          timing.chain_mem_s,
                          timing.chain_xfer_s);
        timing.chain_grow_kernel_s =
            detail_seconds_containing(chain_profiler, "grow_edge_chains_");
        timing.chain_dedup_kernel_s =
            detail_seconds_containing(chain_profiler, "dedup_record_edge_chains");
        timing.num_curvelets = num_curvelets;
        timing.chain_warps_per_block_eff = chain_storage.warp_warps_per_block;
        timing.chain_smem_mode_eff = chain_storage.warp_smem_mode;
        timing.ran_bundle_chain = true;
        chain_profiler.summary();
#if VERBOSE
        std::cout << "Edge chains grown: " << num_curvelets << " curvelets" << std::endl;
#endif

        if (!timing_detail_csv.empty()) {
            chain_profiler.append_detail_csv(timing_detail_csv, timing.run_id);
        }

        std::vector<int> host_chains;
        std::vector<float> host_info;
        if (!gpu_download_compact_chains(chain_storage, host_chains, host_info, num_curvelets)) {
            gpu_curvelet_free_chains(chain_storage);
            gpu_curvelet_free_bundles(bundle_storage);
            gpu_preprocess_free(neighbor_graph);
            return false;
        }

        const unsigned out_w = static_cast<unsigned>(chain_storage.chain_width);
        write_int_array_to_file(out_chain_file, host_chains.data(),
                                static_cast<int>(num_curvelets), static_cast<int>(out_w));

        std::vector<double> host_info_d(host_info.begin(), host_info.end());
        write_double_array_to_file(chain_to_info_filename(out_chain_file), host_info_d.data(),
                                   static_cast<int>(num_curvelets), GPU_CURVELET_INFO_WIDTH);

        gpu_curvelet_free_chains(chain_storage);
        gpu_curvelet_free_bundles(bundle_storage);
    }
    else {
        std::cout << "CSR layout ready: " << neighbor_graph.total_neighbor_pairs
                  << " anchor-neighbor pairs, max number of neighbor edges per anchor = "
                  << neighbor_graph.max_num_of_neighbors << std::endl;
        std::cout << "GPU curve bundle formation requires --neighbor-layout fixed-row." << std::endl;
    }

    gpu_preprocess_free(neighbor_graph);

    timing.wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_t0).count();

    const std::string summary_line = timing_summary_row(timing);
    std::cout << "\nTIMING_CSV " << summary_line << std::endl;
    if (!append_timing_summary_csv(timing_csv, timing)) {
        return false;
    }
    if (!timing_csv.empty()) {
        std::cout << "Wrote timing summary row to " << timing_csv << std::endl;
    }
    if (!timing_detail_csv.empty()) {
        std::cout << "Wrote timing detail rows to " << timing_detail_csv << std::endl;
    }

    (void)curvelet_style;
    (void)out_type;
    return true;
}

int main(int argc, char **argv)
{
    int gpu_id = 0;
    std::string out_file = "chain_gpu.txt";
    std::string timing_csv;
    std::string timing_detail_csv;
    std::string neighbor_degree_csv;
    bool neighbor_degree_only = false;
    CurveletParams params;
    bool show_help = false;

    if (!parse_args(argc, argv, params, out_file, gpu_id, timing_csv, timing_detail_csv,
                    neighbor_degree_csv, neighbor_degree_only, show_help)) {
        return show_help ? 0 : 1;
    }

    const bool ok = run_curvelet_gpu(out_file, gpu_id, params, timing_csv, timing_detail_csv,
                                     neighbor_degree_csv, neighbor_degree_only);
    return ok ? 0 : 1;
}
