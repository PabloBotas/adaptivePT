#ifndef __GPU_RAY_POSITIONING_CUH__
#define __GPU_RAY_POSITIONING_CUH__

#include "patient_parameters.hpp"
#include "special_types.hpp"
#include <vector>

__host__ void virtual_src_to_treatment_plane(const unsigned int& num,
                                             const std::vector<BeamAngles_t>& angles,
                                             const float3& ct_offsets);

__global__ void virtual_src_to_treatment_plane_kernel(const int num,
                                                      const float2 *angles,
                                                      const float3 ct_offsets);
__host__ void treatment_plane_to_virtual_src(Array4<float>& pos,
                                             const Array4<float>& pos2,
                                             const Patient_Parameters_t& pat);
__global__ void treatment_plane_to_virtual_src_kernel(const int num,
                                                      const int nbeams,
                                                      float4* pos_,
                                                      const float4* dir_,
                                                      const float2* angles,
                                                      const float4* plane_dir,
                                                      const float4* plane_pos,
                                                      const float3 ct_offsets);
void correct_offsets(const unsigned int& num,
                     const float3& offsets,
                     const float3& original_offsets);
Array4<float> offset_endpoints(const Array4<float>& orig_endpoints,
                               const float3& offsets,
                               const float3& original_offsets);
__global__ void correct_offsets_kernel(const int num,
                                       const float3 offsets,
                                       const float3 original_offsets);
__global__ void correct_offsets_kernel(const int num,
                                       float4* dev_orig_endpoints,
                                       const float3 offsets,
                                       const float3 original_offsets);

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

// __host__ Planes_t
// get_treatment_planes (const Patient_Parameters_t& pat,
//                       const std::vector<BeamAngles_t>& angles);

__device__ __host__ float3 ext_to_int_coordinates(float3 a);
__device__ __host__ float4 ext_to_int_coordinates(float4 a);
__device__ __host__ float3 int_to_ext_coordinates(float3 a);
__device__ __host__ float4 int_to_ext_coordinates(float4 a);
__device__ float4 rotate(const float4& p, const float& gantry, const float& couch);

#endif
