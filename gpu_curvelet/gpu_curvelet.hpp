#ifndef GPU_CURVELET_HPP
#define GPU_CURVELET_HPP

#include <map>
#include <vector>
#include <deque>
#include <list>
#include <utility> 
#include <iostream>
#include <cmath>

#include "indices.hpp"
#include "curvelet_utils.hpp"
#include "gpu_kernels.hpp"

// cuda api
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

class CurveletGPU
{
protected:
    int _sz_edge_data;
    int _num_edges;
    int _group_max_sz;
    int _max_num_look_edges;

    // > target edge information
    float te_id;
    float te_pt_x;
    float te_pt_y;
    float te_orient;
    float te_grad_mag; 

    // > look edge information
    float le_id;
    float le_pt_x;
    float le_pt_y;
    float le_orient;
    float le_grad_mag;

    float _dx; //< variation of edge in pixel
    float _dt; //< variation of angle in radians
    float _sx; //< sample pixel
    float _st; //< sample angle in radians

    // > number of curves in a bundle in terms of location and orientation variations
    unsigned curves_num_in_bundle_pixel;
    unsigned curves_num_in_bundle_theta;

    // > used for building bundles
    float _max_k;

private:
    //> cpu
    float *bundle_min_ks;
    float *bundle_max_ks;
    float *cmp_bundle_min_ks;
    float *cmp_bundle_max_ks;
    float *intersect_bundle_min_ks;
    float *intersect_bundle_max_ks;

    //> gpu
    int device_id;
    cudaEvent_t start, stop;
    float time_gpu_curvelet;        //> store gpu time

    float *dev_edgeLookList;
    float *dev_bundle_min_ks;
    float *dev_bundle_max_ks;
    float *dev_cmp_bundle_min_ks;
    float *dev_cmp_bundle_max_ks;
    float *dev_intersect_bundle_min_ks;
    float *dev_intersect_bundle_max_ks;
    bool  *dev_hyp_LookEdge; 

    //> RETRIEVE DEBUG
    float *retr_edgeLookList;
    float *dev_retr_edgeLookList;

    bool *hyp_LookEdge;
    unsigned *edge_chain_final;         //< store all final curvelets by edge ids
    unsigned *edge_chain_on_the_fly;    //< store a curvelet candidate on the fly
    unsigned *edge_chain_target;        //< store curvelets w.r.t. target edge

public:

    float *_edgeLookList;
    // > constructor
    CurveletGPU(int device, int &num_edges, int &sz_edge_data, float *edgeLookList, unsigned max_num_look_edges,
                float dx, float dt, float sx, float st, float max_k, unsigned group_max_sz):
                device_id(device), _num_edges(num_edges), _sz_edge_data(sz_edge_data), 
                _dx(dx), _dt(dt), _sx(sx), _st(st), _max_k(max_k), _group_max_sz(group_max_sz), _max_num_look_edges(max_num_look_edges)
    {
        //_edgeLookList = edgeLookList;
        _edgeLookList = new float[(_sz_edge_data+1) * (_max_num_look_edges+1) * _num_edges ];

        //> assign edgeLookList to edgeLookList_fit
        for (unsigned i = 0; i < _num_edges; i++) {
            for (unsigned j = 0; j < (_max_num_look_edges+1)*5; j++) {
                _edgeLookList(i, j) = edgeLookList(i, j);
            }
        }

        //> 1) calculate number of curves in a bundle in terms of pixel and theta
        curves_num_in_bundle_pixel = 2*floor( _dx/_sx + 0.5 ) + 1;
        curves_num_in_bundle_theta = 2*floor( _dt/_st + 0.5 ) + 1;

        //> 2) cpu allocations
        bundle_min_ks = new float[ _max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];
        bundle_max_ks = new float[ _max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];

        cmp_bundle_min_ks = new float[ 2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];
        cmp_bundle_max_ks = new float[ 2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];

        intersect_bundle_min_ks = new float[ curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];
        intersect_bundle_max_ks = new float[ curves_num_in_bundle_pixel * curves_num_in_bundle_theta ];

        hyp_LookEdge = new bool[ _max_num_look_edges ];

        edge_chain_on_the_fly = new unsigned[ _group_max_sz ];
        edge_chain_target     = new unsigned[ (_group_max_sz+1) * _max_num_look_edges ];
        edge_chain_final      = new unsigned[ (_num_edges*_max_num_look_edges) * (_group_max_sz+1) ];

        //> 3) gpu allocations
        cudacheck( cudaMalloc((void**)&dev_edgeLookList,                ((_sz_edge_data+1) * (_max_num_look_edges+1) * _num_edges)*sizeof(float)) );
        cudacheck( cudaMalloc((void**)&dev_bundle_min_ks,               (_max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta) * sizeof(float)) );
        cudacheck( cudaMalloc((void**)&dev_bundle_max_ks,               (_max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta) * sizeof(float)) );
        cudacheck( cudaMalloc((void**)&dev_cmp_bundle_min_ks,           (2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta) * sizeof(float)) );
        cudacheck( cudaMalloc((void**)&dev_cmp_bundle_max_ks,           (2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta) * sizeof(float)) );
        cudacheck( cudaMalloc((void**)&dev_intersect_bundle_min_ks,     (curves_num_in_bundle_pixel * curves_num_in_bundle_theta)*sizeof(float)) );
        cudacheck( cudaMalloc((void**)&dev_intersect_bundle_max_ks,     (curves_num_in_bundle_pixel * curves_num_in_bundle_theta)*sizeof(float)) );
        cudacheck( cudaMalloc((void**)&dev_hyp_LookEdge,                (_max_num_look_edges)*sizeof(bool)) );

        //> RETRIEVE DEBUG
        retr_edgeLookList = new float[(_sz_edge_data+1) * (_max_num_look_edges+1) * _num_edges ];
        cudacheck( cudaMalloc((void**)&dev_retr_edgeLookList, ((_sz_edge_data+1) * (_max_num_look_edges+1) * _num_edges)*sizeof(float)) );

        //> 4) cuda event
        cudacheck( cudaEventCreate(&start) );
        cudacheck( cudaEventCreate(&stop) );
    }

    // > destructor
    ~CurveletGPU();

    void preprocessing();
    void build_curvelets_greedy();
    bool compute_curve_bundle( unsigned te_idx, unsigned le_idx );
    bool bundle_valid_check( unsigned le_idx );
    void retrieve_edge_data_from_edgeLookList(unsigned target_idx, unsigned look_idx, float &id, float &pt_x, float &pt_y, float &orient, float &strength);

    void move_to_cmp_bundle( unsigned cmp_idx, unsigned le_idx, bool rep_by_intersection );
    void bundle_intersection();
    bool bundle_intersection_valid_check();
    bool check_curvelet_exist( unsigned edge_chain_on_the_fly_sz );

    //> for retrieving data from GPU
    void retrieve_data_from_GPU();
};

void CurveletGPU::preprocessing() {
    //> cpu
    for (unsigned i = 0; i < _max_num_look_edges; i++) {
        for (unsigned j = 0; j < (curves_num_in_bundle_pixel * curves_num_in_bundle_theta) ; j++) {
            bundle_min_ks(i, j) = -_max_k;
            bundle_max_ks(i, j) = _max_k;
        }

        hyp_LookEdge[i] = 0;
    }

    //> initialize edge_chain_on_the_fly
    for (unsigned i = 0; i < _group_max_sz; i++) {
        edge_chain_on_the_fly[i] = 0;
    }

    //> initialize edge_chain_target
    for (unsigned i = 0; i < _max_num_look_edges; i++) {
        for (unsigned j = 0; j < (_group_max_sz+1); j++) {
            edge_chain_target(i, j) = 0;
        }
    }

    //> initialize edge_chain_final
    for (unsigned i = 0; i < _num_edges*_max_num_look_edges; i++) {
        for (unsigned j = 0; j < (_group_max_sz+1); j++) {
            edge_chain_final(i, j) = 0;
        }
    }

    //> RETRIEVE DEBUG
    #if 1
    for (unsigned i = 0; i < _num_edges; i++) {
        for (unsigned j = 0; j < (_sz_edge_data+1) * (_max_num_look_edges+1); j++) {
            retr_edgeLookList(i, j) = -1;
        }
    }
    #endif

    //> gpu
    //> 1) assign edgeLookList
    cudacheck( cudaMemcpy(dev_edgeLookList, _edgeLookList, ((_sz_edge_data+1) * (_max_num_look_edges+1) * _num_edges)*sizeof(float), cudaMemcpyHostToDevice) );

    //> 2) initialize other arrays
    cudacheck( cudaMemset(dev_bundle_min_ks,          -_max_k, (_max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta)*sizeof(float)) );
    cudacheck( cudaMemset(dev_bundle_max_ks,           _max_k, (_max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta)*sizeof(float)) );
    cudacheck( cudaMemset(dev_cmp_bundle_min_ks,            0, (2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta)*sizeof(float)) );
    cudacheck( cudaMemset(dev_cmp_bundle_max_ks,            0, (2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta)*sizeof(float)) );
    cudacheck( cudaMemset(dev_intersect_bundle_min_ks,      0, (curves_num_in_bundle_pixel * curves_num_in_bundle_theta)*sizeof(float)) );
    cudacheck( cudaMemset(dev_intersect_bundle_max_ks,      0, (curves_num_in_bundle_pixel * curves_num_in_bundle_theta)*sizeof(float)) );
    cudacheck( cudaMemset(dev_hyp_LookEdge,                 0, (_max_num_look_edges)*sizeof(bool)) );

    std::cout<<"preprocess..."<<std::endl;
    for (unsigned i = 0; i < 3; i++) {
        std::cout<<"Cell ID #"<<retr_edgeLookList(i, 0)<<" Look ID:";
        for (unsigned j = 0; j < _max_num_look_edges; j++) {
            if (retr_edgeLookList(i, (j+1)*5) < 0)
                break;
            std::cout<<retr_edgeLookList(i, (j+1)*5)<<" ";
        }
        std::cout<<std::endl;
    }
}

void CurveletGPU::retrieve_edge_data_from_edgeLookList(unsigned target_idx, unsigned look_idx, float &id, float &pt_x, float &pt_y, float &orient, float &strength)
{
    id       = _edgeLookList(target_idx, look_idx*5);
    pt_x     = _edgeLookList(target_idx, look_idx*5 + 1);
    pt_y     = _edgeLookList(target_idx, look_idx*5 + 2);
    orient   = _edgeLookList(target_idx, look_idx*5 + 3);
    strength = _edgeLookList(target_idx, look_idx*5 + 4);
}

// > highest level main code
void CurveletGPU::build_curvelets_greedy( )
{
    // > Process to build curvelets in a greedy way:
    // > 1) loop over all target edges, and form pairs of curves from target edges and look edges
    // > 2) examine all pairs of curves by intersection to gradually build a curve with maximal grouping size of 7 edges

    //> some variables
    bool valid_bundle_created = false;
    bool valid_bundle_intersection = false;
    int valid_edge_num = 0;
    unsigned edge_chain_target_idx = 0;
    unsigned last_edge_chain_target_idx = 0;
    unsigned edge_chain_lidx = 0;

    unsigned DEBUG_TE_ID = 1;

    //> gpu timing start
	cudacheck( cudaEventRecord(start) );

    //> GPU implemented curvelet building!!
	gpu_build_curvelets_greedy( device_id, _num_edges, _max_num_look_edges,
                                dev_edgeLookList, dev_bundle_min_ks, dev_bundle_max_ks,
                                dev_cmp_bundle_min_ks, dev_cmp_bundle_max_ks, 
                                dev_intersect_bundle_min_ks, dev_intersect_bundle_max_ks,
                                dev_hyp_LookEdge,
                                dev_retr_edgeLookList
                              );
    


    //> gpu timing stop
	cudacheck( cudaEventRecord(stop) );
	cudacheck( cudaEventSynchronize(stop) );
	cudacheck( cudaEventElapsedTime(&time_gpu_curvelet, start, stop) );
    printf(" ## GPU time = %8.4f ms\n", time_gpu_curvelet );

    //> DEBUG!!!!!
    #if 0
    cudacheck( cudaMemcpy(retr_edgeLookList, dev_retr_edgeLookList, ((_sz_edge_data+1) * (_max_num_look_edges+1) * _num_edges)*sizeof(float), cudaMemcpyDeviceToHost) );
    std::cout<<"CPU..."<<std::endl;
    for (unsigned i = 0; i < 3; i++) {
        std::cout<<"Cell ID #"<<_edgeLookList(i, 0)<<" Look ID:";
        for (unsigned j = 0; j < _max_num_look_edges; j++) {
            if (_edgeLookList(i, (j+1)*5) < 0)
                break;
            std::cout<<_edgeLookList(i, (j+1)*5)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<"GPU..."<<std::endl;
    for (unsigned i = 0; i < 3; i++) {
        std::cout<<"Cell ID #"<<retr_edgeLookList(i, 0)<<" Look ID:";
        for (unsigned j = 0; j < _max_num_look_edges; j++) {
            if (retr_edgeLookList(i, (j+1)*5) < 0)
                break;
            std::cout<<retr_edgeLookList(i, (j+1)*5)<<" ";
        }
        std::cout<<std::endl;
    }
    #endif
}

void CurveletGPU::retrieve_data_from_GPU( )
{
    //> retrieve all data from GPU
	cudacheck( cudaMemcpy(bundle_min_ks,               dev_bundle_min_ks,               (_max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta) * sizeof(float),  cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy(bundle_max_ks,               dev_bundle_max_ks,               (_max_num_look_edges * curves_num_in_bundle_pixel * curves_num_in_bundle_theta) * sizeof(float),  cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy(cmp_bundle_min_ks,           dev_cmp_bundle_min_ks,           (2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta) * sizeof(float),   cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy(cmp_bundle_max_ks,           dev_cmp_bundle_max_ks,           (2 * curves_num_in_bundle_pixel * curves_num_in_bundle_theta) * sizeof(float),   cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy(intersect_bundle_min_ks,     dev_intersect_bundle_min_ks,     (curves_num_in_bundle_pixel * curves_num_in_bundle_theta)*sizeof(float),         cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy(intersect_bundle_max_ks,     dev_intersect_bundle_max_ks,     (curves_num_in_bundle_pixel * curves_num_in_bundle_theta)*sizeof(float),         cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy(hyp_LookEdge,                dev_hyp_LookEdge,                (_max_num_look_edges)*sizeof(bool),                                                               cudaMemcpyDeviceToHost) );


    //> print out the retrieved data
    // TODO
}

CurveletGPU::~CurveletGPU() {
    //> cpu
    delete[] _edgeLookList;
    delete[] bundle_min_ks;
    delete[] bundle_max_ks;
    delete[] cmp_bundle_min_ks;
    delete[] cmp_bundle_max_ks;
    delete[] intersect_bundle_min_ks;
    delete[] intersect_bundle_max_ks;
    delete[] hyp_LookEdge;

    delete[] edge_chain_on_the_fly;
    delete[] edge_chain_target;
    delete[] edge_chain_final;

    //> DEBUG
    delete[] retr_edgeLookList;

    //> free memory gpu
    cudacheck( cudaFree(dev_edgeLookList) );
    cudacheck( cudaFree(dev_bundle_min_ks) );
    cudacheck( cudaFree(dev_bundle_max_ks) );
    cudacheck( cudaFree(dev_cmp_bundle_min_ks) );
    cudacheck( cudaFree(dev_cmp_bundle_max_ks) );
    cudacheck( cudaFree(dev_intersect_bundle_min_ks) );
    cudacheck( cudaFree(dev_intersect_bundle_max_ks) );
    cudacheck( cudaFree(dev_hyp_LookEdge) );

    cudacheck( cudaEventDestroy(start) );
    cudacheck( cudaEventDestroy(stop) );
}

#endif // GPU_CURVELET_HPP