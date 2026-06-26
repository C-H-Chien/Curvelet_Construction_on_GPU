#ifndef PARAM_SETTINGS_HPP
#define PARAM_SETTINGS_HPP

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

struct CurveletParams {
    //> The following parameters are the default values for the curvelet construction
    double nrad = 3.5;
    double gap = 1.5;
    double dx = 0.4;
    double dt_deg = 15.0;
    double max_k = 0.3;
    unsigned curvelet_style = 2;
    unsigned group_max_sz = 4;
    unsigned out_type = 0;
    unsigned max_look_edge_num = 0;
    double sx = 0.1;
    double st = 0.08;
    unsigned look_edge_slots = 64;
    std::string edge_file = "eth3d_cables2.txt";
    int edge_data_sz = 4;
};

inline void print_usage(const char *prog)
{
    std::cerr
        << "Usage: " << prog << " [options]\n"
        << "\n"
        << "Options:\n"
        << "  --float, --double          Scalar type (default: double)\n"
        << "  --output <file>            Output chain file (default: chain_cpu.txt)\n"
        << "  --nthreads <N>             OpenMP thread count (default: 1)\n"
        << "  --edge-file <file>         Input edge file in (default: ../test_files/eth3d_cables2.txt)\n"
        << "  --nrad <val>               Neighbor search radius (default: 3.5)\n"
        << "  --dx <val>                 Edge position tolerance in pixels (default: 0.4)\n"
        << "  --dt-deg <val>             Edge angle tolerance in degrees (default: 15)\n"
        << "  --max-k <val>              Maximum curvature (default: 0.3)\n"
        << "  --curvelet-style <N>       Curvelet style (default: 2)\n"
        << "  --group-max-sz <N>         Maximum group size (default: 4)\n"
        << "  --out-type <N>             Output type (default: 0)\n"
        << "  --max-look-edge-num <N>    Initial max look-edge count (default: 0)\n"
        << "  --sx <val>                 Position sampling step (default: 0.1)\n"
        << "  --st <val>                 Angle sampling step in radians (default: 0.08)\n"
        << "  --look-edge-slots <N>      Look-edge table slots per edge (default: 64)\n"
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
                       bool &use_double, std::string &out_file, int &nthreads,
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
        if (std::strcmp(arg, "--float") == 0) {
            use_double = false;
            out_file = "chain_cpu_float.txt";
        }
        else if (std::strcmp(arg, "--double") == 0) {
            use_double = true;
            out_file = "chain_cpu_double.txt";
        }
        else if (std::strcmp(arg, "--output") == 0 && i + 1 < argc) {
            out_file = argv[++i];
        }
        else if (std::strcmp(arg, "--nthreads") == 0 && i + 1 < argc) {
            nthreads = std::atoi(argv[++i]);
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
        else if (std::strcmp(arg, "--max-look-edge-num") == 0 && i + 1 < argc) {
            if (!parse_unsigned_arg(argv[++i], "--max-look-edge-num", params.max_look_edge_num)) return false;
        }
        else if (std::strcmp(arg, "--sx") == 0 && i + 1 < argc) {
            if (!parse_double_arg(argv[++i], "--sx", params.sx)) return false;
        }
        else if (std::strcmp(arg, "--st") == 0 && i + 1 < argc) {
            if (!parse_double_arg(argv[++i], "--st", params.st)) return false;
        }
        else if (std::strcmp(arg, "--look-edge-slots") == 0 && i + 1 < argc) {
            if (!parse_unsigned_arg(argv[++i], "--look-edge-slots", params.look_edge_slots)) return false;
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
