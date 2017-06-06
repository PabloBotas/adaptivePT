#ifndef __GPU_RUN_CUH__
#define __GPU_RUN_CUH__

#include "special_types.hpp"

#include <vector>
#include <string>

void calculateRays(const std::vector<float4>& xbuffer,
                   const std::vector<float4>& vxbuffer,
                   const std::vector<short2>& ixbuffer,
                   const std::vector<BeamAngles_t>& angles,
                   const short* spots_per_beam,
                   const float3& ct_offsets,
                   float4* positions_scorer,
                   float4* directions_scorer,
                   short2* metadata_scorer,
                   float* traces_scorer = NULL);

unsigned int rays_to_device(const std::vector<float4>& xbuffer,
                            const std::vector<float4>& vxbuffer,
                            const std::vector<short2>& ixbuffer,
                            const std::vector<BeamAngles_t>& angles,
                            const float3& ct_offsets);

// __global__ void rays_to_device_kernel(int num, float2 *angles, float3 ct_offsets);

void freeCTMemory();

#endif  // GPMC_CUH
