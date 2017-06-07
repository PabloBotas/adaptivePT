#ifndef __GPU_RAY_POSITIONING_KERNEL_CUH__
#define __GPU_RAY_POSITIONING_KERNEL_CUH__

__global__ void virtual_src_to_treatment_plane_kernel(const int num,
                                                      const float2 *angles,
                                                      const float3 ct_offsets);
__global__ void treatment_plane_to_virtual_src_kernel(const int num,
                                                      const float2 *angles,
                                                      const float3 ct_offsets);
__device__ float4 ray_trace_to_CT_volume(const float4& pos,
                                         const float4& vel);

#endif
