#ifndef __INITIALIZE_RAYS_CUH__
#define __INITIALIZE_RAYS_CUH__

#include "patient_parameters.hpp"

#include <vector>

void tramps_to_virtual_source(const Patient_Parameters_t& pat,
                              std::vector<float4>& xbuffer,
                              std::vector<float4>& vxbuffer,
                              std::vector<short2>& ixbuffer);

float3 iso_to_virtual_src_pos(float z, float2 SAD, float2 spot);

float2 virtual_src_to_iso_pos(float3 p, float2 SAD);

float2 virtual_src_to_iso_pos(float3 pos, float3 cos);

float3 getDirection(float3 pos, float2 spot);

#endif
