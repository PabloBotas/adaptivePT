#ifndef __GPU_RAY_KERNEL_CUH__
#define __GPU_RAY_KERNEL_CUH__

#include "gpu_ray_class.cuh"

__global__ void raytrace_plan_kernel(const short num,
                                     const short *spots_per_beam,
                                     const double4* const orig_endpoints,
                                     double4 *pos_scorer,
                                     float* traces);

__device__ void score_traces(int voxnum, float *traces);

__device__ size_t get_endpoints_index(const short beam_id,
                                      const short spot_id,
                                      const short* spots_per_beam);

#endif
