#ifndef __GPU_INFLUENCE_KERNEL_CUH__
#define __GPU_INFLUENCE_KERNEL_CUH__

#include "gpu_ray_class.cuh"

__global__ void get_influence_kernel(const uint nspots,
                                     const uint nprobes,
                                     float4* influence,
                                     float* spot_weights,
                                     float* inf_volume,
                                     float *new_energies = NULL);

__device__ float wepl_to_point (Ray& ray, float3 stop_point, bool overwrite_energy = false);

#endif
