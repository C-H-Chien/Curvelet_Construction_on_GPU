#include "gpu_kernels.hpp"
#include "indices.hpp"

//> device functions
#include "gpu_dev_utils.cuh"
#include "gpu_dev_curve_bundle.cuh"

template<int _group_max_sz, int _curves_num_in_bundle_pixel, int _curves_num_in_bundle_theta>
__global__
void
gpu_curvelets_kernel(
        int _maxLookEdges, float _max_k, float _dx, float _dt, float _sx, float _st, 
        float* dev_edgeLookList, 
        float* dev_bundle_min_ks, float* dev_bundle_max_ks,
        float* dev_cmp_bundle_min_ks, float* dev_cmp_bundle_max_ks,
        float* dev_intersect_bundle_min_ks, float* dev_intersect_bundle_max_ks,
        bool* dev_hyp_LookEdge,
        float* dev_retr_edgeLookList
)
{
    extern __shared__ double sdata[];

    const int tx = threadIdx.x;
    const int bx = blockIdx.x ;

    //> shared memory ptrs
    float* sEgLookList              = (float*)sdata;
    float* sBundle_min_ks           = sEgLookList                        + (5 * (_maxLookEdges+1));
    float* sBundle_max_ks           = sBundle_min_ks                     + (_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta * _maxLookEdges);
    int* sTarget_edge_chain_lead    = (int*)sBundle_max_ks               + (_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta * _maxLookEdges);
    int* sTarget_edge_chain_nonlead = sTarget_edge_chain_lead            + (_maxLookEdges * _group_max_sz);
    bool* sHyp_LookEdge             = (bool*)(sTarget_edge_chain_nonlead + (_maxLookEdges * _group_max_sz));
    bool* sValid_Edge_Chain_lead    = sHyp_LookEdge                      + (_maxLookEdges);
    bool* sValid_Edge_Chain_nonlead = sValid_Edge_Chain_lead             + (_maxLookEdges);

    //> register for each thread
    /*float te_pt_x, te_pt_y;
    float le_id;
    float le_pt_x;
    float le_pt_y;
    float le_orient;*/

    //> read data from global memory to shared memory
    //> 1) EdgeLookList
    int i = 0;
    #pragma unroll
    for(i = 0; i < 5; i++) {
        sEgLookList[i*(_maxLookEdges) + tx] = dev_edgeLookList[bx * (_maxLookEdges+1)*5 + tx + i*(_maxLookEdges)];
    }
    if(tx < 5) {
        sEgLookList[5*(_maxLookEdges) + tx] = dev_edgeLookList[bx * (_maxLookEdges+1)*5 + tx + 5*(_maxLookEdges)];
    }
    __syncthreads();

    //> 2) bundle_min_ks and bundle_max_ks
    #pragma unroll
    for (unsigned j = 0; j < (_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta) ; j++) {
        sBundle_min_ks(tx, j) = -_max_k;
        sBundle_max_ks(tx, j) = _max_k;
    }

    //> initialization of sTarget_edge_chain_lead
    for (unsigned i = 0; i < _group_max_sz; i++) {
        sTarget_edge_chain_lead(tx, i)    = 0;
        sTarget_edge_chain_nonlead(tx, i) = 0;
    }
    sHyp_LookEdge[tx]             = false;
    sValid_Edge_Chain_lead[tx]    = false;
    sValid_Edge_Chain_nonlead[tx] = false;

    //> declare and allocate registers for each thread
    float rCmp_bundle_min_ks[_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta * 2];
    float rCmp_bundle_max_ks[_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta * 2];
    float rIntersect_bundle_min_ks[_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta];
    float rIntersect_bundle_max_ks[_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta];
    int rEdge_chain_on_the_fly[_group_max_sz];

    //> registered variables
    //bool valid_bundle_created = false;
    bool continue_examine = true;
    bool valid_bundle_intersection = false;
    unsigned edge_chain_lidx = 0;
    unsigned edge_chain_id_match_count = 0;
    float _dir;

    //> initializations
    for (int i = 0; i < _group_max_sz; i++) {
        rEdge_chain_on_the_fly[i] = 0;
    }
    
    //> retrieve target edge id
    //id       = sEgLookList[tx*5];
    //pt_x     = sEgLookList[tx*5 + 1];
    //pt_y     = sEgLookList[tx*5 + 2];
    //orient   = sEgLookList[tx*5 + 3];

    //> CURVELET BUILDING STARTS HERE ...
    #pragma unroll
    for (unsigned f_run = 0; f_run < 2; f_run++) {
        
        //> refresh data
        edge_chain_lidx   = 0;
        sHyp_LookEdge[tx] = false;

        //> calculate the direction from the target edge to the look edge
        _dir = angle_from_pt_to_pt(sEgLookList[1]           /*te_pt_x*/, sEgLookList[2]           /*te_pt_y*/, 
                                   sEgLookList[(tx+1)*5 + 1]/*le_pt_x*/, sEgLookList[(tx+1)*5 + 2]/*le_pt_y*/);

        //> create curve bundles between the target edge and all the look edges
        if ((f_run == 0) && (dot(_dir, sEgLookList[(tx+1)*5 + 3]) > 0)) {
            compute_curve_bundle< _curves_num_in_bundle_pixel, _curves_num_in_bundle_theta >
            ( _dx, _dt, _sx, _st,
            tx,                           //> le_idx 
            sBundle_min_ks, 
            sBundle_max_ks, 
            sHyp_LookEdge,
            sEgLookList[1],               //> te_pt_x
            sEgLookList[2],               //> te_pt_y
            sEgLookList[3],               //> te_orient
            sEgLookList[(tx+1)*5 + 1],    //> le_pt_x
            sEgLookList[(tx+1)*5 + 2],    //> le_pt_y
            sEgLookList[(tx+1)*5 + 3]     //> le_orient
            );
        }
        else if ((f_run == 1) && (dot(_dir, sEgLookList[(tx+1)*5 + 3]) < 0)) {
            compute_curve_bundle< _curves_num_in_bundle_pixel, _curves_num_in_bundle_theta >
            ( _dx, _dt, _sx, _st,
            tx,                           //> le_idx 
            sBundle_min_ks, 
            sBundle_max_ks, 
            sHyp_LookEdge,
            sEgLookList[1],               //> te_pt_x
            sEgLookList[2],               //> te_pt_y
            sEgLookList[3],               //> te_orient
            sEgLookList[(tx+1)*5 + 1],    //> le_pt_x
            sEgLookList[(tx+1)*5 + 2],    //> le_pt_y
            sEgLookList[(tx+1)*5 + 3]     //> le_orient
            );
        }

        //> for each pair-wise curvelet bundle hypothesis formed w.r.t. the target edge,
        //  examine curve bundle intersections with all the rest of the valid look edges
        //> if the edge is valid
        if (sHyp_LookEdge[tx]) {

            //> keep the curvelet bundle
            move_to_cmp_bundle<_curves_num_in_bundle_pixel, _curves_num_in_bundle_theta, 
                            _curves_num_in_bundle_pixel*_curves_num_in_bundle_theta*2, 
                            _curves_num_in_bundle_pixel*_curves_num_in_bundle_theta>
            ( 0, tx, false, rCmp_bundle_min_ks, rCmp_bundle_max_ks, 
            rIntersect_bundle_min_ks, rIntersect_bundle_max_ks,
            sBundle_min_ks, sBundle_max_ks
            );

            //> loop over all the rest of the look edges
            #pragma unroll
            for (int le_remain_idx = 0; le_remain_idx < _maxLookEdges; le_remain_idx++) {
                if (sHyp_LookEdge[le_remain_idx]) {

                    //> refreshing
                    continue_examine = true;

                    //> if the look edge is the edge itself, put the reference edge id to the edge chain list
                    if (tx == le_remain_idx) {

                        //> put the look edge id to the edge chain list
                        rEdge_chain_on_the_fly[ edge_chain_lidx ] = (int)sEgLookList[(tx+1)*5];
                        edge_chain_lidx++;

                        //> continue looping if the edge chain is not filled entirely by edge ids
                        if (edge_chain_lidx == (_group_max_sz-1)) {
                            break;
                        }
                        else {
                            continue_examine = false;
                        }
                    }

                    //> continue to merge curve bundles with other look edges
                    if (continue_examine) {

                        //> keep the curvelet bundle
                        move_to_cmp_bundle<_curves_num_in_bundle_pixel, _curves_num_in_bundle_theta, 
                                        _curves_num_in_bundle_pixel*_curves_num_in_bundle_theta*2, 
                                        _curves_num_in_bundle_pixel*_curves_num_in_bundle_theta>
                        ( 1, le_remain_idx, false, rCmp_bundle_min_ks, rCmp_bundle_max_ks, 
                        rIntersect_bundle_min_ks, rIntersect_bundle_max_ks,
                        sBundle_min_ks, sBundle_max_ks
                        );

                        //> find bundle intersections: the absolute min and max curvatures
                        bundle_intersection<_curves_num_in_bundle_pixel, _curves_num_in_bundle_theta,
                                            _curves_num_in_bundle_pixel*_curves_num_in_bundle_theta*2, 
                                            _curves_num_in_bundle_pixel*_curves_num_in_bundle_theta>
                        ( rCmp_bundle_min_ks, rCmp_bundle_max_ks, rIntersect_bundle_min_ks, rIntersect_bundle_max_ks );

                        //> check whether the intersected bundle is valid or not
                        valid_bundle_intersection = bundle_intersection_valid_check<_curves_num_in_bundle_pixel, _curves_num_in_bundle_theta,
                                                                                    _curves_num_in_bundle_pixel*_curves_num_in_bundle_theta>
                                                    ( rIntersect_bundle_min_ks, rIntersect_bundle_max_ks );

                        //> if the intersection is valid
                        if (valid_bundle_intersection) {

                            //> put the look edge id to the edge chain list
                            rEdge_chain_on_the_fly[ edge_chain_lidx ] = (int)sEgLookList[(le_remain_idx+1)*5];
                            edge_chain_lidx++;

                            //> change the cmp_bundle set by putting in the intersected bundle
                            //> keep the curvelet bundle
                            move_to_cmp_bundle<_curves_num_in_bundle_pixel, _curves_num_in_bundle_theta, 
                                            _curves_num_in_bundle_pixel*_curves_num_in_bundle_theta*2, 
                                            _curves_num_in_bundle_pixel*_curves_num_in_bundle_theta>
                            ( 0, le_remain_idx, true, rCmp_bundle_min_ks, rCmp_bundle_max_ks, 
                            rIntersect_bundle_min_ks, rIntersect_bundle_max_ks,
                            sBundle_min_ks, sBundle_max_ks
                            );

                        }
                    }

                    //> continue looping if the edge chain is not filled entirely by edge ids
                    if (edge_chain_lidx == (_group_max_sz-1)) {
                        break;
                    }

                }
            }        
        }

        //> refresh
        //> reuse the memory from sHyp_LookEdge
        sHyp_LookEdge[tx] = false;

        //> push every edge_chain_on_the_fly stored in thread registers to sTarget_edge_chain_lead stored in shared memory
        if (edge_chain_lidx > 2) {
            #pragma unroll
            for (unsigned i = 0; i < _group_max_sz; i++) {
                if (f_run == 0) {
                    sTarget_edge_chain_lead(tx, i) = rEdge_chain_on_the_fly[i];
                }
                else if (f_run == 1) {
                    sTarget_edge_chain_nonlead(tx, i) = rEdge_chain_on_the_fly[i];
                }
            }

            //> set true if the local edge_chain_on_the_fly is a candidate
            sHyp_LookEdge[tx] = true;
        }
    
        //> for each thread, compare self rEdge_chain_on_the_fly with others by accessing the sTarget_edge_chain_lead
        bool uniqueness = true;
        if (sHyp_LookEdge[tx]) {
            #pragma unroll
            for (unsigned le = 0; le < _maxLookEdges; le++) {

                //> if the edge_chain_on_the_fly is a candidate in the target edge chain pool
                if (sHyp_LookEdge[le]) {

                    //> reset the match edgel id counter
                    edge_chain_id_match_count = 0;

                    //> loop over all edgel ids
                    for (unsigned j = 0; j < _group_max_sz; j++) {

                        //> compare curvelet id one by one
                        if ((f_run == 0) && (rEdge_chain_on_the_fly[j] == sTarget_edge_chain_lead(le, j))) {
                            edge_chain_id_match_count++;
                        }
                        else if ((f_run == 1) && (rEdge_chain_on_the_fly[j] == sTarget_edge_chain_nonlead(le, j))) {
                            edge_chain_id_match_count++;
                        }
                    }
                    if ((edge_chain_id_match_count == _group_max_sz) && (uniqueness) && (tx != le)) {
                        if (tx < le) {
                            if (f_run == 0) {
                                sValid_Edge_Chain_lead[tx] = true;
                            }
                            else if (f_run == 1) {
                                sValid_Edge_Chain_nonlead[tx] = true;
                            }
                            
                        }
                        else {
                            uniqueness = false;
                        }     
                    }
                }
            }
        }

    }
}

void gpu_build_curvelets_greedy(
        int device_id,
        int num_edges, int maxLookEdges,  
        float* dev_edgeLookList, float* dev_bundle_min_ks, float* dev_bundle_max_ks,
        float* dev_cmp_bundle_min_ks, float* dev_cmp_bundle_max_ks,
        float* dev_intersect_bundle_min_ks, float* dev_intersect_bundle_max_ks,
        bool* dev_hyp_LookEdge,
        float* dev_retr_edgeLookList
)
{
    //> kernel parameters
    float    PI = 3.14159265358979323846;
    float    _max_k = 0.3;
    float    _dx = 0.4;
    float    _dt = (15.0/180.0)*PI;
    float    _sx = 0.1;
    float    _st = 0.08;
    const int _group_max_sz = 7;
    const int _curves_num_in_bundle_pixel = 9; //2*floor( _dx/_sx + 0.5 ) + 1;
    const int _curves_num_in_bundle_theta = 7; //2*floor( _dt/_st + 0.5 ) + 1;

    const int _num_edges = num_edges;
    int _maxLookEdges = maxLookEdges;

    //> intermidiate constants used as template variables
    const float _dx2 = _dx * _dx;

    //> configurations
    const int batchCount = _num_edges;
    const int thread_x = _maxLookEdges;
    dim3 grid(batchCount, 1, 1);
    dim3 threads(thread_x, 1, 1);

    //> declare shared memory storage
    int shmem = 0;
    shmem += sizeof(float)    * (5 * (_maxLookEdges+1));                                                       //> edgeLookList for each warp
    shmem += sizeof(float)    * (_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta * _maxLookEdges);   //> dev_bundle_min_ks
    shmem += sizeof(float)    * (_curves_num_in_bundle_pixel * _curves_num_in_bundle_theta * _maxLookEdges);   //> dev_bundle_max_ks
    shmem += sizeof(int)      * (_maxLookEdges * (_group_max_sz));                                             //> target edge chain pool leading
    shmem += sizeof(int)      * (_maxLookEdges * (_group_max_sz));                                             //> target edge chain pool nonleading
    shmem += sizeof(bool)     * (_maxLookEdges);                                                               //> dev_hyp_LookEdge
    shmem += sizeof(bool)     * (_maxLookEdges);                                                               //> sValidEdgeChain leading
    shmem += sizeof(bool)     * (_maxLookEdges);                                                               //> sValidEdgeChain non leanding

    // get max. dynamic shared memory on the GPU
    int nthreads_max, shmem_max = 0;
    cudacheck( cudaDeviceGetAttribute(&nthreads_max, cudaDevAttrMaxThreadsPerBlock, device_id) );
    #if CUDA_VERSION >= 9000
    cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlockOptin, device_id) );
    if (shmem <= shmem_max) {
        cudacheck( cudaFuncSetAttribute(gpu_curvelets_kernel<_group_max_sz, _curves_num_in_bundle_pixel, _curves_num_in_bundle_theta>, cudaFuncAttributeMaxDynamicSharedMemorySize, shmem) );
    }
    #else
    cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlock, device_id) );
    #endif    // CUDA_VERSION >= 9000

    if ( shmem > shmem_max ) {
        printf("error: kernel %s requires too many threads or too much shared memory\n", __func__);
    }

    void *kernel_args[] = {&_maxLookEdges, &_max_k, &_dx, &_dt, &_sx, &_st, &dev_edgeLookList, &dev_bundle_min_ks, &dev_bundle_max_ks, &dev_cmp_bundle_min_ks, &dev_cmp_bundle_max_ks, 
                           &dev_intersect_bundle_min_ks, &dev_intersect_bundle_max_ks, &dev_hyp_LookEdge, &dev_retr_edgeLookList};

    cudacheck( cudaLaunchKernel((void*)gpu_curvelets_kernel<_group_max_sz, _curves_num_in_bundle_pixel, _curves_num_in_bundle_theta>, grid, threads, kernel_args, shmem, NULL) );
}
