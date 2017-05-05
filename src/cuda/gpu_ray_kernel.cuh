#ifndef __GPU_RAY_KERNEL_CUH__
#define __GPU_RAY_KERNEL_CUH__

__global__ void calculateRays_kernel(const int num,
                                     float4 *scorer,
                                     float* traces,
                                     const short *spots_per_beam);

__device__ unsigned int get_endpoints_index(const short beam_id,
                                            const short spot_id,
                                            const short* spots_per_beam);

#endif
