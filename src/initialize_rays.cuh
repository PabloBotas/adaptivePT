#ifndef __INITIALIZE_RAYS_CUH__
#define __INITIALIZE_RAYS_CUH__

#include "patient_parameters.hpp"
#include "helper_math.h"

#include <vector>

void create_virtual_source_buffers(const Patient_Parameters_t& pat,
                              std::vector<float4>& xbuffer,
                              std::vector<float4>& vxbuffer,
                              std::vector<short2>& ixbuffer);

void set_treatment_plane_buffers_ct_space(const Patient_Parameters_t& pat,
                                          const Array4<float>& endpoints,
                                          const Array4<float>& init_pos,
                                          std::vector<float4>& xbuffer,
                                          std::vector<float4>& vxbuffer,
                                          std::vector<short2>& ixbuffer);

float3 iso_to_virtual_src_pos(float z, float2 SAD, float2 spot);

float2 virtual_src_to_iso_pos(float3 p, float2 SAD);
float2 virtual_src_to_iso_pos(float3 pos, float3 cos);
void virtual_src_to_iso_pos(Array4<float>& pos, SAD_t SAD);

float3 getDirection(float3 pos, float2 spot);

short2 get_beam_spot_id (size_t num, const std::vector<short>& spots_per_field);
#endif
