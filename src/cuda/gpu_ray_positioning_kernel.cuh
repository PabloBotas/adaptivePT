#ifndef __GPU_RAY_POSITIONING_KERNEL_CUH__
#define __GPU_RAY_POSITIONING_KERNEL_CUH__

__global__ void rays_to_delivery_plane(const int num,
                                       const float2 *angles,
                                       const float3 ct_offsets);

__device__ float4 ray_trace_to_CT_volume(const float4& pos,
                                         const float4& vel);

#endif
