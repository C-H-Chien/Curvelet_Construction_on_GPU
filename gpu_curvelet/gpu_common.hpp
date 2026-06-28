#ifndef GPU_COMMON_HPP
#define GPU_COMMON_HPP

#include <cuda_runtime.h>
#include <cstdio>

#define cudacheck(a) do { \
    cudaError_t e = (a); \
    if (e != cudaSuccess) { \
        fprintf(stderr, "\033[1;31mError in %s:%d %s\033[0m\n", \
                __func__, __LINE__, cudaGetErrorString(e)); \
    } \
} while (0)

inline int div_up(int a, int b)
{
    return (a + b - 1) / b;
}

#endif // GPU_COMMON_HPP
