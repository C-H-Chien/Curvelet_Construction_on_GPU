#include <iostream>
#include <cmath>
#include <math.h>
#include <string>
#include <vector>

#include "data_io.hpp"
#include "param_settings.hpp"
#include "preprocess.hpp"
#include "cpu_curvelet.hpp"
#include "timer.hpp"

template<typename T>
const char* scalar_name()
{
    return (sizeof(T) == sizeof(double)) ? "double" : "float";
}

template<typename T>
bool run_curvelet(const std::string &out_chain_file, int nthreads, const CurveletParams &params)
{
    const T nrad = T(params.nrad);
    const T dx = T(params.dx);
    const T dt = T(params.dt_deg / 180.0) * T(M_PI);
    const T max_k = T(params.max_k);
    const unsigned curvelet_style = params.curvelet_style;
    const unsigned group_max_sz = params.group_max_sz;
    const unsigned out_type = params.out_type;
    unsigned max_LookEdgeNum = params.max_look_edge_num;
    const T sx = T(params.sx);
    const T st = T(params.st);
    const unsigned LOOK_EDGE_SLOTS = params.look_edge_slots;

    const std::string &edge_file = params.edge_file;
    int edge_data_sz = params.edge_data_sz;

    std::cout << "Using scalar type: " << scalar_name<T>() << std::endl;

    StepTimer timer;
    timer.start();

    std::vector<T> TOED_edges;
    if (!read_TO_edges_from_file(edge_file, edge_data_sz, TOED_edges)) {
        return false;
    }
    int edge_num = static_cast<int>(TOED_edges.size() / edge_data_sz);
    timer.lap("read TO edges");

    // ------------------------------------------------------------------
    //> preprocessings ...
    const int edgeLookList_src_stride = (edge_data_sz + 1) * (int)LOOK_EDGE_SLOTS;
    T *edgeLookList = new T[edgeLookList_src_stride * edge_num];

    //> 1) construct edge maps and edge lists
    edgeNeighborList<T> edgeLookListObj( edge_num, edge_data_sz, TOED_edges.data(), group_max_sz, nrad, LOOK_EDGE_SLOTS );
    timer.lap("construct edge map");

    edgeLookListObj.init_edgeLookList( edgeLookList );
    timer.lap("init edgeLookList");

    edgeLookListObj.create_edgeLookList( edgeLookList, max_LookEdgeNum );
    timer.lap("create edgeLookList");

    //> 2) build curvelet
    CurveletCPU<T> CurveletCPU_obj( edge_num, edge_data_sz, edgeLookList, max_LookEdgeNum,
                                    edgeLookList_src_stride, dx, dt, sx, st, max_k, group_max_sz, nthreads,
                                    TOED_edges.data(), nrad, T(1));
    timer.lap("CurveletCPU setup (copy edgeLookList)");

    CurveletCPU_obj.build_curvelets_greedy();
    timer.lap("build_curvelets_greedy");

    const unsigned out_h = CurveletCPU_obj.num_curvelets();
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
    timer.lap("write output");
    timer.summary();

    delete[] edgeLookList;
    delete[] out_chain;
    delete[] out_info;

    (void)curvelet_style;
    (void)out_type;
    return true;
}

int main(int argc, char **argv)
{
    int nthreads = 1;
    bool use_double = true;
    std::string out_file = "chain_cpu.txt";
    CurveletParams params;
    bool show_help = false;

    if (!parse_args(argc, argv, params, use_double, out_file, nthreads, show_help)) {
        return show_help ? 0 : 1;
    }

    const bool ok = use_double
        ? run_curvelet<double>(out_file, nthreads, params)
        : run_curvelet<float>(out_file, nthreads, params);

    return ok ? 0 : 1;
}
