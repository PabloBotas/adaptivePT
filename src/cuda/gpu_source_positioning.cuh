#ifndef __GPU_RAY_POSITIONING_CUH__
#define __GPU_RAY_POSITIONING_CUH__

#include "special_types.hpp"
#include <vector>

__host__ void virtual_src_to_treatment_plane(const unsigned int& num,
                                             const std::vector<BeamAngles_t>& angles,
                                             const float3& ct_offsets);

__global__ void virtual_src_to_treatment_plane_kernel(const int num,
                                                      const float2 *angles,
                                                      const float3 ct_offsets);
__global__ void treatment_plane_to_virtual_src_kernel(const int num,
                                                      const float2 *angles,
                                                      const float3 ct_offsets);
__device__ float4 ray_trace_to_CT_volume(const float4& p,
                                         const float4& v);
__device__ __host__ float4 ray_trace_to_CT_volume(const float4& p,
                                                  const float4& v,
                                                  const int3 nvox,
                                                  const float3 dvox);
__device__ __host__ float3 ray_trace_to_CT_volume(const float3& p,
                                                  const float3& v,
                                                  const int3 nvox,
                                                  const float3 dvox);

__host__ Array4<float>
get_treatment_planes (const std::vector<BeamAngles_t>& angles);

__device__ __host__ float3 ext_to_int_coordinates(float3 a);
__device__ __host__ float4 ext_to_int_coordinates(float4 a);
__device__ __host__ float3 int_to_ext_coordinates(float3 a);
__device__ __host__ float4 int_to_ext_coordinates(float4 a);
__device__ float4 rotate(const float4& p, const float& gantry, const float& couch);
__host__ Vector4_t<float> rotate(const Vector4_t<float>& p, const float& gantry, const float& couch);

#endif
