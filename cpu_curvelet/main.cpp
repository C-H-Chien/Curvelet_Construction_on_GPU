#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>

#include "data_io.hpp"
#include "param_settings.hpp"
#include "preprocess.hpp"
#include "cpu_curvelet.hpp"
#include "timer.hpp"

namespace {

struct RunTimingSummary {
    std::string run_id;
    std::string edge_file;
    std::string scalar_type;
    int nthreads = 1;
    int num_edges = 0;
    unsigned max_candidates = 0;

    unsigned max_num_of_neighbors = 0;
    int total_neighbor_pairs = 0;
    double mean_neighbors = 0.0;
    unsigned valid_pairs = 0;
    unsigned num_curvelets = 0;

    double preprocess_s = 0.0;
    double bundle_bench_s = 0.0;       //> standalone form_pairwise_bundles
    double curvelet_build_s = 0.0;     //> build_curvelets_greedy wall
    double curvelet_pairwise_s = 0.0;  //> dir filter + bundle transport (max over threads)
    double curvelet_dir_filter_s = 0.0;
    double curvelet_bundle_transport_s = 0.0;
    double curvelet_grow_s = 0.0;
    double curvelet_dedup_s = 0.0;
    double wall_s = 0.0;
};

const char *timing_summary_header()
{
    return
        "run_id,edge_file,scalar_type,nthreads,num_edges,max_candidates,"
        "max_num_of_neighbors,total_neighbor_pairs,mean_neighbors,"
        "valid_pairs,num_curvelets,"
        "preprocess_s,bundle_bench_s,curvelet_build_s,"
        "curvelet_pairwise_s,curvelet_dir_filter_s,curvelet_bundle_transport_s,"
        "curvelet_grow_s,curvelet_dedup_s,wall_s";
}

std::string timing_summary_row(const RunTimingSummary &r)
{
    std::ostringstream oss;
    oss << csv_escape(r.run_id) << ','
        << csv_escape(r.edge_file) << ','
        << csv_escape(r.scalar_type) << ','
        << r.nthreads << ','
        << r.num_edges << ','
        << r.max_candidates << ','
        << r.max_num_of_neighbors << ','
        << r.total_neighbor_pairs << ','
        << csv_format_seconds(r.mean_neighbors) << ','
        << r.valid_pairs << ','
        << r.num_curvelets << ','
        << csv_format_seconds(r.preprocess_s) << ','
        << csv_format_seconds(r.bundle_bench_s) << ','
        << csv_format_seconds(r.curvelet_build_s) << ','
        << csv_format_seconds(r.curvelet_pairwise_s) << ','
        << csv_format_seconds(r.curvelet_dir_filter_s) << ','
        << csv_format_seconds(r.curvelet_bundle_transport_s) << ','
        << csv_format_seconds(r.curvelet_grow_s) << ','
        << csv_format_seconds(r.curvelet_dedup_s) << ','
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
    std::string base = edge_file;
    const auto slash = base.find_last_of("/\\");
    if (slash != std::string::npos) {
        base = base.substr(slash + 1);
    }
    return base + "_" + std::to_string(ms);
}

} // namespace

template<typename T>
const char* scalar_name()
{
    return (sizeof(T) == sizeof(double)) ? "double" : "float";
}

template<typename T>
bool run_curvelet(const std::string &out_chain_file, int nthreads, const CurveletParams &params,
                  const std::string &timing_csv, const std::string &timing_detail_csv)
{
    const T nrad = T(params.nrad);
    const T dx = T(params.dx);
    const T dt = T(params.dt_deg / 180.0) * T(M_PI);
    const T max_k = T(params.max_k);
    const unsigned curvelet_style = params.curvelet_style;
    const unsigned group_max_sz = params.group_max_sz;
    const unsigned out_type = params.out_type;
    const T sx = T(params.sx);
    const T st = T(params.st);

    const std::string &edge_file = params.edge_file;
    int edge_data_sz = params.edge_data_sz;

    // std::cout << "Using scalar type: " << scalar_name<T>() << std::endl;

    std::vector<T> TOED_edges;
    if (!read_TO_edges_from_file(edge_file, edge_data_sz, TOED_edges)) {
        return false;
    }
    int edge_num = static_cast<int>(TOED_edges.size() / edge_data_sz);

    RunTimingSummary timing;
    timing.run_id = make_run_id(edge_file);
    timing.edge_file = edge_file;
    timing.scalar_type = scalar_name<T>();
    timing.nthreads = nthreads;
    timing.num_edges = edge_num;
    timing.max_candidates = params.max_candidates;

    const auto wall_t0 = std::chrono::steady_clock::now();

    StepTimer timer;
    timer.start();

    //> Preprocess: identify neighbor edges for each anchor edge and structure the neighbor graph in compressed sparse rows (CSR) form
    CPUNeighborGraph<T> csr_graph;
    if (!build_neighbor_csr_graph( edge_num, edge_data_sz, TOED_edges.data(), nrad, 3u, params.max_candidates, nthreads, csr_graph)) {
        return false;
    }
    timing.preprocess_s = timer.lap("build CSR neighbor graph");
    timing.max_num_of_neighbors = csr_graph.max_num_of_neighbors;
    timing.total_neighbor_pairs = csr_graph.total_neighbor_pairs;
    timing.mean_neighbors = (edge_num > 0)
        ? static_cast<double>(csr_graph.total_neighbor_pairs) / static_cast<double>(edge_num)
        : 0.0;

    //> Curve bundle formation: pairwise transport (standalone benchmark; also redone inside build)
    CurveletCPU<T> CurveletCPU_obj( edge_num, edge_data_sz, csr_graph, dx, dt, sx, st, max_k, group_max_sz, nthreads, TOED_edges.data(), nrad );

    CPUCurveletFormationResult bundle_result;
    if (!CurveletCPU_obj.form_pairwise_bundles(bundle_result)) {
        return false;
    }
    timing.bundle_bench_s = timer.lap("form_pairwise_bundles");
    timing.valid_pairs = bundle_result.valid_pairs;
    std::cout << "Pairwise curve bundles formed: " << bundle_result.valid_pairs << " valid pairs" << std::endl;

    //> Chain growth + dedup/record (internal printout splits phase 1 vs phase 2)
    CPUCurveletBuildTiming build_timing;
    CurveletCPU_obj.build_curvelets_greedy(&build_timing);
    timing.curvelet_build_s = timer.lap("build_curvelets_greedy (pairwise+growth+dedup; see phase breakdown above)");
    timing.curvelet_dir_filter_s = build_timing.direction_filter_s;
    timing.curvelet_bundle_transport_s = build_timing.bundle_transport_s;
    timing.curvelet_pairwise_s = build_timing.direction_filter_s + build_timing.bundle_transport_s;
    timing.curvelet_grow_s = build_timing.chain_growth_s;
    timing.curvelet_dedup_s = build_timing.dedup_s;
    //> Prefer omp wall from inside the parallel region when available.
    if (build_timing.wall_s > 0.0) {
        timing.curvelet_build_s = build_timing.wall_s;
    }

    timer.add("curvelet_dir_filter (max over threads)", build_timing.direction_filter_s);
    timer.add("curvelet_bundle_transport (max over threads)", build_timing.bundle_transport_s);
    timer.add("curvelet_grow (max over threads)", build_timing.chain_growth_s);
    timer.add("curvelet_dedup (max over threads)", build_timing.dedup_s);
    timer.summary();

    const unsigned out_h = CurveletCPU_obj.num_curvelets();
    timing.num_curvelets = out_h;
    const unsigned out_w = CurveletCPU_obj.chain_width();
    const unsigned info_w = 10;
    std::cout<<"(out_h, out_w) = ("<<out_h<<", "<<out_w<<")"<<std::endl;

    int *out_chain = new int[out_h * out_w];
    for (unsigned i = 0; i < out_h; i++) {
        for (unsigned j = 0; j < out_w; j++) {
            out_chain[i * out_w + j] = (int)CurveletCPU_obj._edge_chain_final[i * out_w + j];
        }
    }
    write_int_array_to_file(out_chain_file, out_chain, out_h, out_w);

    double *out_info = new double[info_w * out_h];
    const unsigned info_stride = CurveletCPU_obj._max_curvelets;
    for (unsigned row = 0; row < out_h; row++) {
        for (unsigned col = 0; col < info_w; col++) {
            out_info[row * info_w + col] = CurveletCPU_obj._curvelet_info[col * info_stride + row];
        }
    }
    write_double_array_to_file(chain_to_info_filename(out_chain_file), out_info, out_h, info_w);

    delete[] out_chain;
    delete[] out_info;

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
        if (!timer.append_detail_csv(timing_detail_csv, timing.run_id)) {
            return false;
        }
        std::cout << "Wrote timing detail rows to " << timing_detail_csv << std::endl;
    }

    (void)curvelet_style;
    (void)out_type;
    return true;
}

int main(int argc, char **argv)
{
    int nthreads = 1;
    bool use_double = true;
    std::string out_file = "chain_cpu.txt";
    std::string timing_csv;
    std::string timing_detail_csv;
    CurveletParams params;
    bool show_help = false;

    if (!parse_args(argc, argv, params, use_double, out_file, nthreads,
                    timing_csv, timing_detail_csv, show_help)) {
        return show_help ? 0 : 1;
    }

    const bool ok = use_double
        ? run_curvelet<double>(out_file, nthreads, params, timing_csv, timing_detail_csv)
        : run_curvelet<float>(out_file, nthreads, params, timing_csv, timing_detail_csv);

    return ok ? 0 : 1;
}
