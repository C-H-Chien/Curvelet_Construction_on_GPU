#ifndef PARAM_SETTINGS_HPP
#define PARAM_SETTINGS_HPP

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

struct CurveletParams {
    double nrad = 3.5;
    double dx = 0.4;
    double dt_deg = 15.0;
    double max_k = 0.3;
    unsigned curvelet_style = 2;
    unsigned group_max_sz = 4;
    unsigned out_type = 0;
    double sx = 0.1;
    double st = 0.08;
    std::string edge_file = "eth3d_cables2.txt";
    int edge_data_sz = 4;
    std::string csr_strategy = "two-pass";
    std::string csr_discover_mode = "thread";
    std::string neighbor_layout = "fixed-row";
    unsigned max_candidates = 32;
    int neighbor_count_threads = 1;
    int neighbor_fill_threads = 1;
    int neighbor_stage_threads = 1;
    int neighbor_warps_per_block = 1;
    std::string fixed_row_build = "warp";
    int bundle_warps_per_block = 1; //> pairwise curve-bundle formation: warps/block
    int chain_warps_per_block = 2; //> TODO: tune this parameter
    //> Shared-memory cache for warp growth: auto | none | lane
    //> auto  = prefer lane workspace in shared (may reduce warps/block)
    //> none  = no shared cache (mode 0); keep requested warps/block
    //> lane  = per-lane growth workspace in shared (mode 1); keep requested warps/block
    std::string chain_smem_mode = "auto";
    int dedup_threads_per_block = 128;

    //> Angle tolerance in radians (kernels use radians; CLI stores degrees).
    float dt_rad() const { return static_cast<float>(dt_deg * M_PI / 180.0); }
};

inline void print_usage(const char *prog)
{
    std::cerr
        << "Usage: " << prog << " [options]\n"
        << "\n"
        << "Options:\n"
        << "  --output <file>            Output chain file (default: chain_gpu.txt)\n"
        << "  --device <N>               CUDA device id (default: 0)\n"
        << "  --edge-file <file>         Input edge file in test_files\n"
        << "  --nrad <val>               Neighbor search radius (default: 3.5)\n"
        << "  --dx <val>                 Edge position tolerance in pixels (default: 0.4)\n"
        << "  --dt-deg <val>             Edge angle tolerance in degrees (default: 15)\n"
        << "  --max-k <val>              Maximum curvature (default: 0.3)\n"
        << "  --curvelet-style <N>       Curvelet style (default: 2)\n"
        << "  --group-max-sz <N>         Maximum group size (default: 4)\n"
        << "  --out-type <N>             Output type (default: 0)\n"
        << "  --sx <val>                 Position sampling step (default: 0.1)\n"
        << "  --st <val>                 Angle sampling step in radians (default: 0.08)\n"
        << "  --csr-strategy <mode>      CSR build: single-pass | two-pass (default: two-pass)\n"
        << "  --csr-discover <mode>      CSR discover, using warp or thread per anchor edge: thread | warp (default: thread)\n"
        << "  --neighbor-layout <mode>   Neighbor storage: csr | fixed-row (default: fixed-row)\n"
        << "  --fixed-row-build <mode>   Fixed-row build: warp | stage (default: warp)\n"
        << "  --neighbor-warps-per-block <N>  Warp-per-anchor discover: warps/block (default: 1)\n"
        << "  --bundle-warps-per-block <N>  Pairwise bundle formation: warps/block (default: 4)\n"
        << "  --chain-warps-per-block <N>  Warp-per-anchor chain growth: warps/block (default: 4)\n"
        << "  --chain-smem-mode <mode>   Warp shared cache: auto | none | lane (default: auto)\n"
        << "  --dedup-threads-per-block <N>  Deduplication threads/block (default: 128)\n"
        << "  --max-candidates <N>       Max neighbors staged per anchor (default: 64)\n"
        << "  --neighbor-count-threads <N>  Two-pass count kernel threads/block (default: 1)\n"
        << "  --neighbor-fill-threads <N>   Two-pass fill kernel threads/block (default: 1)\n"
        << "  --neighbor-stage-threads <N>  Single-pass stage/compact threads/block (default: 1)\n"
        << "  --edge-data-sz <N>         Values per edge in input file (default: 4)\n"
        << "  --help                     Show this help message\n";
}

inline bool parse_unsigned_arg(const char *arg, const char *name, unsigned &value)
{
    char *end = nullptr;
    const unsigned long parsed = std::strtoul(arg, &end, 10);
    if (end == arg || *end != '\0') {
        std::cerr << "Error: invalid value for " << name << ": " << arg << std::endl;
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

inline bool parse_int_arg(const char *arg, const char *name, int &value)
{
    char *end = nullptr;
    const long parsed = std::strtol(arg, &end, 10);
    if (end == arg || *end != '\0') {
        std::cerr << "Error: invalid value for " << name << ": " << arg << std::endl;
        return false;
    }
    value = static_cast<int>(parsed);
    return true;
}

inline bool parse_double_arg(const char *arg, const char *name, double &value)
{
    char *end = nullptr;
    const double parsed = std::strtod(arg, &end);
    if (end == arg || *end != '\0') {
        std::cerr << "Error: invalid value for " << name << ": " << arg << std::endl;
        return false;
    }
    value = parsed;
    return true;
}

inline bool parse_args(int argc, char **argv, CurveletParams &params,
                       std::string &out_file, int &gpu_id,
                       bool &show_help)
{
    show_help = false;
    for (int i = 1; i < argc; i++) {
        const char *arg = argv[i];
        if (std::strcmp(arg, "--help") == 0 || std::strcmp(arg, "-h") == 0) {
            print_usage(argv[0]);
            show_help = true;
            return false;
        }
        else if (std::strcmp(arg, "--output") == 0 && i + 1 < argc) {
            out_file = argv[++i];
        }
        else if (std::strcmp(arg, "--device") == 0 && i + 1 < argc) {
            if (!parse_int_arg(argv[++i], "--device", gpu_id)) return false;
        }
        else if ((std::strcmp(arg, "--input") == 0 || std::strcmp(arg, "--edge-file") == 0) && i + 1 < argc) {
            params.edge_file = argv[++i];
        }
        else if (std::strcmp(arg, "--nrad") == 0 && i + 1 < argc) {
            if (!parse_double_arg(argv[++i], "--nrad", params.nrad)) return false;
        }
        else if (std::strcmp(arg, "--dx") == 0 && i + 1 < argc) {
            if (!parse_double_arg(argv[++i], "--dx", params.dx)) return false;
        }
        else if (std::strcmp(arg, "--dt-deg") == 0 && i + 1 < argc) {
            if (!parse_double_arg(argv[++i], "--dt-deg", params.dt_deg)) return false;
        }
        else if (std::strcmp(arg, "--max-k") == 0 && i + 1 < argc) {
            if (!parse_double_arg(argv[++i], "--max-k", params.max_k)) return false;
        }
        else if (std::strcmp(arg, "--curvelet-style") == 0 && i + 1 < argc) {
            if (!parse_unsigned_arg(argv[++i], "--curvelet-style", params.curvelet_style)) return false;
        }
        else if (std::strcmp(arg, "--group-max-sz") == 0 && i + 1 < argc) {
            if (!parse_unsigned_arg(argv[++i], "--group-max-sz", params.group_max_sz)) return false;
        }
        else if (std::strcmp(arg, "--out-type") == 0 && i + 1 < argc) {
            if (!parse_unsigned_arg(argv[++i], "--out-type", params.out_type)) return false;
        }
        else if (std::strcmp(arg, "--sx") == 0 && i + 1 < argc) {
            if (!parse_double_arg(argv[++i], "--sx", params.sx)) return false;
        }
        else if (std::strcmp(arg, "--st") == 0 && i + 1 < argc) {
            if (!parse_double_arg(argv[++i], "--st", params.st)) return false;
        }
        else if (std::strcmp(arg, "--csr-strategy") == 0 && i + 1 < argc) {
            params.csr_strategy = argv[++i];
        }
        else if (std::strcmp(arg, "--csr-discover") == 0 && i + 1 < argc) {
            params.csr_discover_mode = argv[++i];
        }
        else if (std::strcmp(arg, "--neighbor-layout") == 0 && i + 1 < argc) {
            params.neighbor_layout = argv[++i];
        }
        else if (std::strcmp(arg, "--fixed-row-build") == 0 && i + 1 < argc) {
            params.fixed_row_build = argv[++i];
        }
        else if (std::strcmp(arg, "--neighbor-warps-per-block") == 0 && i + 1 < argc) {
            if (!parse_int_arg(argv[++i], "--neighbor-warps-per-block", params.neighbor_warps_per_block)) return false;
        }
        else if (std::strcmp(arg, "--bundle-warps-per-block") == 0 && i + 1 < argc) {
            if (!parse_int_arg(argv[++i], "--bundle-warps-per-block", params.bundle_warps_per_block)) return false;
        }
        else if (std::strcmp(arg, "--chain-warps-per-block") == 0 && i + 1 < argc) {
            if (!parse_int_arg(argv[++i], "--chain-warps-per-block", params.chain_warps_per_block)) return false;
        }
        else if (std::strcmp(arg, "--chain-smem-mode") == 0 && i + 1 < argc) {
            params.chain_smem_mode = argv[++i];
        }
        else if (std::strcmp(arg, "--dedup-threads-per-block") == 0 && i + 1 < argc) {
            if (!parse_int_arg(argv[++i], "--dedup-threads-per-block", params.dedup_threads_per_block)) return false;
        }
        else if (std::strcmp(arg, "--max-candidates") == 0 && i + 1 < argc) {
            if (!parse_unsigned_arg(argv[++i], "--max-candidates", params.max_candidates)) return false;
        }
        else if (std::strcmp(arg, "--neighbor-count-threads") == 0 && i + 1 < argc) {
            if (!parse_int_arg(argv[++i], "--neighbor-count-threads", params.neighbor_count_threads)) return false;
        }
        else if (std::strcmp(arg, "--neighbor-fill-threads") == 0 && i + 1 < argc) {
            if (!parse_int_arg(argv[++i], "--neighbor-fill-threads", params.neighbor_fill_threads)) return false;
        }
        else if (std::strcmp(arg, "--neighbor-stage-threads") == 0 && i + 1 < argc) {
            if (!parse_int_arg(argv[++i], "--neighbor-stage-threads", params.neighbor_stage_threads)) return false;
        }
        else if (std::strcmp(arg, "--edge-data-sz") == 0 && i + 1 < argc) {
            if (!parse_int_arg(argv[++i], "--edge-data-sz", params.edge_data_sz)) return false;
        }
        else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            print_usage(argv[0]);
            return false;
        }
    }
    return true;
}

#endif
