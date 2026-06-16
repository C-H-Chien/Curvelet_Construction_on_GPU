
/*****************************************************************************
// file: main
// author: Chiang-Heng Chien
// date: 06/14/2022
//       An algorithm to form curvelet from input third order edge list
******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <ctime>

#include "preprocess.hpp"
#include "cpu_curvelet.hpp"

/*************************************************************
Usage: 
 [chain, info] = form_curvelet_mex(edgeinfo, nrows, ncols,...
 rad, gap, dx, dt, token_len, max_k, cvlet_type, ...
 max_size_to_goup, output_type)
Input:
        edgeinfo: nx4 array storing the position, orientation
 and magnitude information;
        nrows: height of the image;
        ncols: width of the image;
        rad: radius of the grouping neighborhood around each
 edgel;
        gap: distance between two consecutive edges in a link;
        dx: position uncertainty at the reference edgel;
        dt: orientation uncertainty at the reference edgel;
        token_len:
        max_k: maximum curvature of curves in the curve bundle;
        cvlet_type: if 0, form regular anchor centered curvelet;
 if 1, anchor centered bidirectional; if 2, anchor leading bidirectional;
 if 3, ENO style (anchor leading or trailing but in the same direction).
        max_size_to_goup: the maximum numbers of edges to group;
        output_type: if 0, out put curvelet map. if 1, out put
 the curve fragment map. if 2, out put the poly arc map.
 *************************************************************/

template<typename T>
void read_TO_edges_from_file(std::string filename, T *rd_data, int first_dim, int second_dim)
{
    std::cout<<"reading data from a file "<<filename<<std::endl;
    std::string in_file_name = "../test_files/";
    in_file_name.append(filename);
    std::fstream in_file;
    T data;
    int j = 0, i = 0;

    in_file.open(in_file_name, std::ios_base::in);
    if (!in_file) {
        std::cerr << "input read file not existed!\n";
    }
    else {
        while (in_file >> data) {
            rd_data[i * second_dim + j] = data;
            j++;
            if (j == second_dim) {
                j = 0;
                i++;
            }
        }
    }
}

void write_int_array_to_file(std::string filename, int *wr_data, int first_dim, int second_dim)
{
    std::cout<<"writing data to a file "<<filename<<" ..."<<std::endl;
    std::string out_file_name = "../outputs/";
    out_file_name.append(filename);
	std::ofstream out_file;
    out_file.open(out_file_name);
    if ( !out_file.is_open() )
      std::cout<<"write data file cannot be opened!"<<std::endl;

	for (int i = 0; i < first_dim; i++) {
		for (int j = 0; j < second_dim; j++) {
            out_file << wr_data[i * second_dim + j] <<"\t";
		}
		out_file << "\n";
	}

    out_file.close();
}

template<typename T>
const char* scalar_name()
{
    return (sizeof(T) == sizeof(double)) ? "double" : "float";
}

template<typename T>
void run_curvelet(const std::string &out_chain_file)
{
    const T PI = T(3.14159265358979323846);

    // -- input settings (match original_code / cpu_curvelet_omp for comparison)
    int height = 800;
    int width = 800;
    T nrad = T(3.5);
    T gap = T(1.5);
    T dx = T(0.4);
    T dt = T(15.0/180.0)*PI;
    T token_len = T(1);
    T max_k = T(0.3);
    unsigned curvelet_style = 2;
    unsigned group_max_sz = 4;
    unsigned out_type = 0;
    unsigned max_LookEdgeNum = 0;
    T sx = T(0.1);
    T st = T(0.08);
    const unsigned LOOK_EDGE_SLOTS = 32;

    int edge_num = 3965;
    int edge_data_sz = 4;

    std::cout << "Using scalar type: " << scalar_name<T>() << std::endl;

    T *TOED_edges;
    TOED_edges = new T[edge_num * edge_data_sz];

    // read third-order edges from file
    read_TO_edges_from_file("TO_edges_ABC_0006_thresh1.txt", TOED_edges, edge_num, edge_data_sz);

    //std::cout<<TOED_edges[0]<<"  "<<TOED_edges[1]<<"  "<<TOED_edges[2]<<std::endl;

    // ------------------------------------------------------------------
    // > preprocessings ...
    const unsigned edge_look_stride = (edge_data_sz + 1) * LOOK_EDGE_SLOTS;
    T *edgeLookList = new T[edge_look_stride * edge_num];

    //edgeNeighborList( int &neighbor_sz, const T &rad):

    // > 1) construct edge maps and edge lists
    edgeNeighborList<T> edgeLookListObj( width, height, edge_num, edge_data_sz, TOED_edges, group_max_sz, nrad, LOOK_EDGE_SLOTS );
    edgeLookListObj.init_edgeLookList( edgeLookList );
    edgeLookListObj.create_edgeLookList( edgeLookList, max_LookEdgeNum );

    //std::cout<<"maxmial number of look edges (from main) = "<<max_LookEdgeNum<<std::endl;

    // > 2) build curvelet
    CurveletCPU<T> CurveletCPU_obj( edge_num, edge_data_sz, edgeLookList, max_LookEdgeNum, dx, dt, sx, st, max_k, group_max_sz);
    CurveletCPU_obj.preprocessing();
    CurveletCPU_obj.build_curvelets_greedy();

    const unsigned out_h = CurveletCPU_obj.num_curvelets();
    const unsigned out_w = CurveletCPU_obj.chain_width();
    std::cout<<"(out_h, out_w) = ("<<out_h<<", "<<out_w<<")"<<std::endl;

    int *out_chain = new int[out_h * out_w];
    for (unsigned i = 0; i < out_h; i++) {
        for (unsigned j = 0; j < out_w; j++) {
            out_chain[i * out_w + j] = (int)CurveletCPU_obj._edge_chain_final[i * out_w + j];
        }
    }
    write_int_array_to_file(out_chain_file, out_chain, out_h, out_w);

    delete[] TOED_edges;
    delete[] edgeLookList;
    delete[] out_chain;

    (void)gap;
    (void)token_len;
    (void)curvelet_style;
    (void)out_type;
}

int main(int argc, char **argv)
{
    bool use_double = true;
    std::string out_file = "chain_cpu.txt";

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--float") == 0) {
            use_double = false;
            out_file = "chain_cpu_float.txt";
        } else if (std::strcmp(argv[i], "--double") == 0) {
            use_double = true;
            out_file = "chain_cpu_double.txt";
        } else if (std::strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            out_file = argv[++i];
        }
    }

    if (use_double)
        run_curvelet<double>(out_file);
    else
        run_curvelet<float>(out_file);

    return 0;
}
