#ifndef __INITIALIZE_RAYS_CUH__
#define __INITIALIZE_RAYS_CUH__

#include "patient_parameters.hpp"

#include <vector>

void init_rays(const Patient_Parameters_t& pat,
               std::vector<float4>& xbuffer,
               std::vector<float4>& vxbuffer,
               std::vector<short2>& ixbuffer);

float3 adjust_to_internal_coordinates(float3 a);

float3 getTanslatedPosition(float z, float2 SAD, float2 spot);

float3 getDirection(float3 pos, float2 spot);

#endif
