#ifndef __GPU_INFLUENCE_KERNEL_CUH__
#define __GPU_INFLUENCE_KERNEL_CUH__

#include "gpu_ray_class.cuh"

__global__ void get_influence_kernel(const short num,
                                     const double4* const endpos,
                                     double4* influence,
                                     float* spot_weights,
                                     float* inf_volume);

__device__ double wepl_to_point (Ray& ray, double3 stop_point, bool overwrite_energy = false);

__device__ double get_influence (double const energy, double const wepl_r, double const wepl_d);

#endif
