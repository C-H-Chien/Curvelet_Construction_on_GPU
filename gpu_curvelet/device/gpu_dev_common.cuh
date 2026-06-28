#ifndef GPU_DEV_COMMON_CUH
#define GPU_DEV_COMMON_CUH

#include <cuda_runtime.h>

constexpr int kEdgeFields = 4;

__device__ __forceinline__ long long pack_cell_key(int x, int y)
{
    //> pack the cell key as a 64-bit integer
    return (static_cast<long long>(y) << 32) | (static_cast<unsigned int>(x));
}

__device__ __forceinline__ float sq_dist_dev(float x1, float y1, float x2, float y2)
{
    const float dx = x1 - x2;
    const float dy = y1 - y2;
    return dx * dx + dy * dy;
}

//> Use binary search to discover the neighbor edges for each anchor edge
__device__ __forceinline__ int lower_bound_cell_keys(const long long *keys, int n, long long key)
{
    int lo = 0;
    int hi = n;
    while (lo < hi) {
        const int mid = lo + ((hi - lo) >> 1);
        if (keys[mid] < key) {
            lo = mid + 1;
        } 
        else {
            hi = mid;
        }
    }
    return lo;
}

#endif // GPU_DEV_COMMON_CUH
