
/*****************************************************************************
// file: form_curvelet_mex.cxx
// author: Xiaoyan Li
// date: 01/19/2015
//       An algorithm to form curvelet from input edgemap
******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <fstream>
#include <string.h>
#include <vector>
#include <ctime>
#include "Array.hpp"
#include "form_curvelet_process.hpp"
#include "curvelet_utils.hpp"

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

void read_TO_edges_from_file(std::string filename, double *rd_data, int first_dim, int second_dim)
{
#define rd_data(i, j) rd_data[(i) * second_dim + (j)]
    std::cout<<"reading data from a file "<<filename<<std::endl;
    std::string in_file_name = "../test_files/";
    in_file_name.append(filename);
    std::fstream in_file;
    double data;
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

void write_double_array_to_file(std::string filename, double *wr_data, int first_dim, int second_dim)
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
    const double PI = 3.14159265358979323846;

    // -- settings
    //int height = 481; 
    //int width = 321; 
    int height = 480; 
    int width = 752; 
    double nrad = 3.5;
    double gap = 1.5;
    double dx = 0.4;
    double dt = (15.0/180.0)*PI; 
    double token_len = 1;
    double max_k = 0.3;
    unsigned curvelet_style = 3;
    unsigned group_max_sz = 7;
    unsigned out_type = 0;

    int edge_num = 72043;
    //int edge_num = 100;
    int edge_data_sz = 4;
    double *TOED_edges;
    TOED_edges = new double[edge_num * edge_data_sz];

    // read third-order edges from file
    //read_TO_edges_from_file("TO_edges.txt", TOED_edges, edge_num, edge_data_sz);
    read_TO_edges_from_file("TO_edges_digit1_T.txt", TOED_edges, edge_num, edge_data_sz);
    //read_TO_edges_from_file("TO_edges_EuRoC.txt", TOED_edges, edge_num, edge_data_sz);

    // construct and assign the subpixel edge list
    arrayd edgeinfo; 
    edgeinfo._data = TOED_edges;
    
    int h = edge_num; 
    edgeinfo.set_h(h);
    int w = edge_data_sz; 
    edgeinfo.set_w(w);

    //std::cout<<edgeinfo._data[100]<<"\t"<<edgeinfo._data[101]<<std::endl;
    //std::cout<<edgeinfo._data[28793]<<"\t"<<edgeinfo._data[28794]<<std::endl;

    unsigned output_type = out_type;
    
    // assign initial settings
    unsigned cvlet_type = curvelet_style;
    unsigned max_size_to_group = group_max_sz;
    bool bCentered_grouping = cvlet_type==0 || cvlet_type==1, bBidirectional_grouping = cvlet_type==0 || cvlet_type==2;
    
    form_curvelet_process curvelet_pro(edgeinfo, unsigned(height),unsigned(width),
                                       nrad, gap,
                                       dx, dt,
                                       token_len, max_k,
                                       max_size_to_group,
                                       bCentered_grouping, bBidirectional_grouping);
    printf("process is constructed\n");
    
    
    curvelet_pro.execute();
    printf("process is executed\n");
    
    // create memory for output
    unsigned out_h,out_w, info_w;
    
    curvelet_pro.get_output_size(out_h,out_w,output_type);

    std::cout<<"(out_h, out_w) = ("<<out_h<<", "<<out_w<<")"<<std::endl;
    int *out_chain = new int[out_h * out_w];

    if(output_type==0)
        info_w = 10;
    else if(output_type==1)
        info_w = 1;
    else
        info_w = 12;
    //const mwSize ds2[2] = {mwSize(outh), mwSize(infow)};
    //const unsigned ds2[2] = {out_h, info_w};

    double *out_info = new double[out_h * info_w];
    arrayi chain;
    //pl[0] = mxCreateNumericArray(2,ds1,mxINT32_CLASS, mxREAL);
    //chain._data = (int*) mxGetData(pl[0]);
    chain._data = out_chain;
    chain.set_h(out_h);
    chain.set_w(out_w);
    arrayd info;
    info._data = out_info;

    info.set_h(out_h);
    info.set_w(info_w);
    
    curvelet_pro.get_output_arrary( chain, info, output_type );
/*
    for (int j = 0; j < 10; j++) {
        for (int i = 0; i < out_w; i++) {
            std::cout<< std::right << std::setw(3) << chain._data[out_h*i+j] << ", ";
        }
        std::cout<<std::endl;
    }
*/
    //write_int_array_to_file("chain.txt", out_chain, out_h, out_w);
    //write_int_array_to_file("chain.txt", chain._data, out_h, out_w);
    //write_double_array_to_file("info.txt", info._data, info_w, out_h);

    delete[] TOED_edges;
    delete[] out_info;
    delete[] out_chain;

}

