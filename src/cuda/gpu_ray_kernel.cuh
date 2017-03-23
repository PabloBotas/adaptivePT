#ifndef __GPU_RAY_KERNEL_CUH__
#define __GPU_RAY_KERNEL_CUH__

__global__ void calculateRays_kernel(int num, float *scorer);
__global__ void rays_to_device_kernel(int num, float2 angles, float3 ct_offsets);

#endif
