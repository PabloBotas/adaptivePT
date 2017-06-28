#include "initialize_rays.cuh"

#include "tramp.hpp"
#include "gpu_source_positioning.cuh"
#include "spot.hpp"
#include "helper_math.h"

#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <cuda.h>

#define MeV2eV 1e6

void create_virtual_source_buffers(const Patient_Parameters_t& pat,
                              std::vector<float4>& xbuffer,
                              std::vector<float4>& vxbuffer,
                              std::vector<short2>& ixbuffer)
//  initialize particle buffer
{
    std::cout << "Gantry: " << pat.machine << std::endl;
    for (unsigned int ibeam = 0; ibeam < pat.beam_names.size(); ibeam++)
    {
        std::cout << "Loading beam: " << pat.beam_names.at(ibeam) << std::endl;
        std::cout << "    Beam #:     " << ibeam << std::endl;
        std::cout << "    Tramp file: " << pat.tramp_files.at(ibeam) << std::endl;

        Tramp_t src(pat.tramp_files.at(ibeam), pat.machine);
        src.z = - pat.isocenter_to_beam_distance.at(ibeam); // cm

        if(pat.range_shifters[ibeam].exists)
        {
            RangeShifter_Dims_t rs = pat.range_shifters[ibeam];
            std::cout << "    Range shifter thickness:    " << rs.thick << " cm" << std::endl;
            for (size_t i = 0; i < src.wepls.size(); i++)
            {
                src.wepls.at(i) -= rs.wepl;
            }
        }
        std::cout << "    Source Z plane: " << src.z << " cm" << std::endl;

        float2 SAD = make_float2(pat.virtualSAD.a, pat.virtualSAD.b); // cm

        xbuffer.reserve(src.nspots);
        vxbuffer.reserve(src.nspots);
        ixbuffer.reserve(src.nspots);
        for(unsigned int i=0; i < src.nspots; i++)
        { // LOOP OVER SPOTS
            Spot_t& spot = src.spots[i];

            float3 pos  = iso_to_virtual_src_pos(src.z, SAD, make_float2(spot.x, spot.y)); // cm
            float3 dCos = getDirection(pos, make_float2(spot.x, spot.y));
            float energy = src.energies_internal.at(i)*MeV2eV; // eV
            float wepl   = src.wepls.at(i);                    // cm

            xbuffer.push_back( make_float4(pos.x, pos.y, pos.z, wepl) );
            vxbuffer.push_back( make_float4(dCos.x, dCos.y, dCos.z, energy) );
            ixbuffer.push_back( make_short2(ibeam, i) );
        }
    }
}

void create_treatment_plane_buffers (const Patient_Parameters_t& pat,
                                     const Array4<float>& endpoints,
                                     const Array4<float>& init_pos,
                                     std::vector<float4>& xbuffer,
                                     std::vector<float4>& vxbuffer,
                                     std::vector<short2>& ixbuffer)
{
    size_t s = endpoints.size();
    xbuffer.resize(s);
    vxbuffer.resize(s);
    ixbuffer.resize(s);
    for (size_t i = 0; i < s; i++)
    {
        float3 start = make_float3(init_pos.at(i).x, init_pos.at(i).y, init_pos.at(i).z);
        float3 end   = make_float3(endpoints.at(i).x, endpoints.at(i).y, endpoints.at(i).z);
        float3 dir   = end - start;
        float3 dCos  = dir/length(dir);
        float wepl   = init_pos.at(i).w;
        float energy = endpoints.at(i).w;
        short2 meta  = get_beam_spot_id(i, pat.spots_per_field);

        int3 nvox   = make_int3(pat.ct.n.x, pat.ct.n.y, pat.ct.n.z);
        float3 dvox = make_float3(pat.ct.d.x, pat.ct.d.y, pat.ct.d.z);
        float3 start2 = ray_trace_to_CT_volume(start, dCos, nvox, dvox);

        xbuffer.at(i)  = make_float4(start2, wepl);
        vxbuffer.at(i) = make_float4(dCos, energy);
        ixbuffer.at(i) = meta;

        // printf("%d - 0 - %f %f %f - %f %f %f - %f %f %f - %f %f %f - %f %f %f\n", 
        //        i, xbuffer[i].x, xbuffer[i].y, xbuffer[i].z,
        //        vxbuffer[i].x, vxbuffer[i].y, vxbuffer[i].z,
        //        xbuffer[i].x, xbuffer[i].y, xbuffer[i].z,
        //        vxbuffer[i].x, vxbuffer[i].y, vxbuffer[i].z,
        //        end.x, end.y, end.z);
    }
}

float3 iso_to_virtual_src_pos(float z, float2 SAD, float2 spot)
{
    float3 p;
    p.x = ((SAD.x + z) / SAD.x) * spot.x;
    p.y = ((SAD.y + z) / SAD.y) * spot.y;
    p.z = z;
    return p;
}

float2 virtual_src_to_iso_pos(float3 pos, float2 SAD)
{
    float2 spot;
    spot.x = pos.x * SAD.x / (SAD.x + pos.z);
    spot.y = pos.y * SAD.y / (SAD.y + pos.z);
    return spot;
}

void virtual_src_to_iso_pos(Array4<float>& pos, SAD_t SAD)
{
    for (size_t i = 0; i < pos.size(); i++)
    {
        pos.at(i).x = pos.at(i).x * SAD.a / (SAD.a + pos.at(i).z);
        pos.at(i).y = pos.at(i).y * SAD.b / (SAD.b + pos.at(i).z);
    }
}

float2 virtual_src_to_iso_pos(float3 pos, float3 cos)
{
    float2 spot;
    spot.x = pos.x + abs(pos.z)*cos.x / sqrt(1-cos.x*cos.x);
    spot.y = pos.y + abs(pos.z)*cos.y / sqrt(1-cos.y*cos.y);
    return spot;
}

float3 getDirection(float3 pos, float2 spot)
{
    float3 dCos;
    float a = (spot.x-pos.x)/abs(pos.z);
    float b = (spot.y-pos.y)/abs(pos.z);
    float norm = sqrt(a*a + b*b + 1.f);
    dCos.x = a/norm;
    dCos.y = b/norm;

    float temp = 1.0f - dCos.x*dCos.x - dCos.y*dCos.y;
    if(temp < 0)
    {
        std::cerr << "Something went wrong calculating direction cosines:\n";
        std::cerr << "    Pos  x:  " << pos.x  << "\n";
        std::cerr << "    Pos  y:  " << pos.y  << "\n";
        std::cerr << "    Pos  z:  " << pos.z  << "\n";
        std::cerr << "    spot x:  " << spot.x << "\n";
        std::cerr << "    spot y:  " << spot.y << std::endl;
        exit(EXIT_FAILURE);
    };
    dCos.z = sqrt(temp);
    return dCos;
}

short2 get_beam_spot_id (size_t num, const std::vector<short>& spots_per_field)
{
    size_t beamid = 0;
    for (; beamid < spots_per_field.size(); beamid++)
    {
        if (num >= (size_t)spots_per_field.at(beamid))
            num -= (size_t)spots_per_field.at(beamid);
        else
            break;
    }
    size_t spotid = num;
    return make_short2(beamid, spotid);
}
