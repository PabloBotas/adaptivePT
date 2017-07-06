#ifndef __GPU_RAY_KERNEL_CUH__
#define __GPU_RAY_KERNEL_CUH__

#include "gpu_ray_class.cuh"
#include <curand_kernel.h>

template<class T>
__global__ void raytrace_plan_kernel(const short num,
                                     const short *spots_per_beam,
                                     const float4* const orig_endpoints,
                                     float4 *pos_scorer,
                                     T* traces);

template<class T>
__device__ void score_traces(int voxnum, curandState& localState, const short& beamid, const short& spotid, T *traces);

__device__ size_t get_endpoints_index(const short beam_id,
                                      const short spot_id,
                                      const short* spots_per_beam);

#endif
