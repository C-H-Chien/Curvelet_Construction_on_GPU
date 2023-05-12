#ifndef GPU_KERNELS_HPP
#define GPU_KERNELS_HPP

#include<stdio.h>
#include<assert.h>
#include<cuda.h>
#include<cuda_runtime_api.h>

// cuda error check
#define cudacheck( a )  do { \
                            cudaError_t e = a; \
                            if(e != cudaSuccess) { \
                                printf("\033[1;31m"); \
                                printf("Error in %s:%d %s\n", __func__, __LINE__, cudaGetErrorString(e)); \
                                printf("\033[0m"); \
                            }\
                        } while(0)

//#ifdef __cplusplus
//extern "C" {
//#endif

//> single precision curvelet building
void gpu_build_curvelets_greedy(
        int device_id,
        int num_edges, int maxLookEdges,  
        float* dev_edgeLookList, float* dev_bundle_min_ks, float* dev_bundle_max_ks,
        float* dev_cmp_bundle_min_ks, float* dev_cmp_bundle_max_ks,
        float* dev_intersect_bundle_min_ks, float* dev_intersect_bundle_max_ks,
        bool* dev_hyp_LookEdge,
        float* dev_retr_edgeLookList
);

//#ifdef __cplusplus
//    }
//#endif

#endif // GPU_KERNELS_HPP
