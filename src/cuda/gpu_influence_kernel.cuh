#ifndef __GPU_INFLUENCE_KERNEL_CUH__
#define __GPU_INFLUENCE_KERNEL_CUH__

#include "gpu_ray_class.cuh"

__global__ void get_influence_kernel(const short num,
                                     const short* spots_per_field,
                                     const double4* const endpos,
                                     double4* influence);

__device__ double wepl_to_point (Ray& ray, double3 stop_point);

__device__ double get_influence (double wepl_r, double wepl_d);

#endif
