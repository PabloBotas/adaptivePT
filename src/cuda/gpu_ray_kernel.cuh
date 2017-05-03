#ifndef __GPU_RAY_KERNEL_CUH__
#define __GPU_RAY_KERNEL_CUH__

__global__ void calculateRays_kernel(const int num,
                                     float4 *scorer,
                                     const short *spots_per_beam);

__global__ void rays_to_device_kernel(const int num,
                                      float2 *angles,
                                      const float3 ct_offsets);

__device__ unsigned int find_scorer_index(const short beam_id,
                                          const short spot_id,
                                          const short* spots_per_beam);

#endif
