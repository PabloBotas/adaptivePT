#ifndef __INITIALIZE_RAYS_CUH__
#define __INITIALIZE_RAYS_CUH__

#include "patient_parameters.hpp"

#include <vector>

void create_virtual_source_buffers(const Patient_Parameters_t& pat,
                              std::vector<float4>& xbuffer,
                              std::vector<float4>& vxbuffer,
                              std::vector<short2>& ixbuffer);

void create_treatment_plane_buffers(const Patient_Parameters_t& pat,
                                    const std::vector< Vector4_t<float> >& endpoints,
                                    const std::vector< Vector4_t<float> >& init_pos,
                                    std::vector<float4>& xbuffer,
                                    std::vector<float4>& vxbuffer,
                                    std::vector<short2>& ixbuffer);

float3 iso_to_virtual_src_pos(float z, float2 SAD, float2 spot);

float2 virtual_src_to_iso_pos(float3 p, float2 SAD);
float2 virtual_src_to_iso_pos(float3 pos, float3 cos);

float3 getDirection(float3 pos, float2 spot);
float3 getDirection(float3 dir);

short2 get_beam_spot_id (size_t num, const std::vector<short> spots_per_field);
#endif
