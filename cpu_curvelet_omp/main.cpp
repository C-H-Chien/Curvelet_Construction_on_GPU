
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
#include <string.h>
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

void read_TO_edges_from_file(std::string filename, float *rd_data, int first_dim, int second_dim)
{
#define rd_data(i, j) rd_data[(i) * second_dim + (j)]
    std::cout<<"reading data from a file "<<filename<<std::endl;
    std::string in_file_name = "../test_files/";
    in_file_name.append(filename);
    std::fstream in_file;
    float data;
    int j = 0, i = 0;

    in_file.open(in_file_name, std::ios_base::in);
    if (!in_file) {
        std::cerr << "input read file not existed!\n";
    }
    else {
        while (in_file >> data) {
            rd_data(i, j) = data;
            j++;
            if (j == second_dim) {
                j = 0;
                i++;
            }
        }
    }
#undef rd_data
}

void write_double_array_to_file(std::string filename, float *wr_data, int first_dim, int second_dim)
{
#define wr_data(i, j) wr_data[(i) * second_dim + (j)]

    std::cout<<"writing data to a file "<<filename<<" ..."<<std::endl;
    std::string out_file_name = "../test_files/";
    out_file_name.append(filename);
	std::ofstream out_file;
    out_file.open(out_file_name);
    if ( !out_file.is_open() )
      std::cout<<"write data file cannot be opened!"<<std::endl;

	for (int i = 0; i < first_dim; i++) {
		for (int j = 0; j < second_dim; j++) {
			out_file << wr_data(i, j) <<"\t";
		}
		out_file << "\n";
	}

    out_file.close();
#undef wr_data
}

void write_int_array_to_file(std::string filename, int *wr_data, int first_dim, int second_dim)
{
#define wr_data(i, j) wr_data[(i) * second_dim + (j)]

    std::cout<<"writing data to a file "<<filename<<" ..."<<std::endl;
    std::string out_file_name = "../test_files/";
    out_file_name.append(filename);
	std::ofstream out_file;
    out_file.open(out_file_name);
    if ( !out_file.is_open() )
      std::cout<<"write data file cannot be opened!"<<std::endl;

	for (int i = 0; i < first_dim; i++) {
		for (int j = 0; j < second_dim; j++) {
			//out_file << wr_data(i, j) <<"\t";
            out_file << wr_data[i + first_dim * j] <<"\t";
		}
		out_file << "\n";
	}

    out_file.close();
#undef wr_data
}

int main(int argc, char **argv)
{
    int nthreads = 8;
    

    const float PI = 3.14159265358979323846;

    // -- input settings
    int height = 481; 
    int width = 321; 
    //int height = 28; 
    //int width = 28; 
    float nrad = 3.5;
    float gap = 1.5;
    float dx = 0.4;
    float dt = (15.0/180.0)*PI; 
    float token_len = 1;
    float max_k = 0.3;
    unsigned curvelet_style = 3;
    unsigned group_max_sz = 7;
    unsigned out_type = 0;
    unsigned max_LookEdgeNum = 0;
    float sx = 0.1;
    float st = 0.08;

    //int edge_num = 28793;
    int edge_num = 31042;
    //int edge_num = 100;
    int edge_data_sz = 4;
    float *TOED_edges;
    TOED_edges = new float[edge_num * edge_data_sz];

    // read third-order edges from file
    read_TO_edges_from_file("TO_edges.txt", TOED_edges, edge_num, edge_data_sz);
    //read_TO_edges_from_file("TO_edges_digit1.txt", TOED_edges, edge_num, edge_data_sz);

    //std::cout<<TOED_edges[0]<<"  "<<TOED_edges[1]<<"  "<<TOED_edges[2]<<std::endl;

    // ------------------------------------------------------------------
    // > preprocessings ...
    float *edgeLookList = new float[(edge_data_sz+1) * 32 * edge_num ];

    //edgeNeighborList( int &neighbor_sz, const float &rad):

    // > 1) construct edge maps and edge lists
    edgeNeighborList edgeLookListObj( width, height, edge_num, edge_data_sz, TOED_edges, group_max_sz, nrad );
    edgeLookListObj.init_edgeLookList( edgeLookList );
    edgeLookListObj.create_edgeLookList( edgeLookList, max_LookEdgeNum );

    //std::cout<<"maxmial number of look edges (from main) = "<<max_LookEdgeNum<<std::endl;

    // > 2) build curvelet. Now using float only.
    CurveletCPU<float> CurveletCPU_fp32( edge_num, edge_data_sz, edgeLookList, max_LookEdgeNum, dx, dt, sx, st, max_k, group_max_sz, nthreads);
    CurveletCPU_fp32.build_curvelets_greedy();

    //edgeMap _edgeMap(width, height, edge_num, edge_data_sz, TOED_edges);
    //_edgeMap.print_map();


    delete[] TOED_edges;
    delete[] edgeLookList;
}

