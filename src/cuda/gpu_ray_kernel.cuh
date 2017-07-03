#ifndef __GPU_RAY_KERNEL_CUH__
#define __GPU_RAY_KERNEL_CUH__

template<class T>
__global__ void raytrace_plan_kernel(const short num,
                                     const short *spots_per_beam,
                                     const float4* const orig_endpoints,
                                     float4 *pos_scorer,
                                     T* traces);

template<class T>
__device__ score_traces(const int thread, Ray ray, T *traces, int voxnum);

__device__ size_t get_endpoints_index(const short beam_id,
                                      const short spot_id,
                                      const short* spots_per_beam);

#endif
