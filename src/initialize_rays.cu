#include "initialize_rays.cuh"
#include "tramp.hpp"
#include "spot.hpp"

#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <cuda.h>

#define MeV2eV 1e6

void tramps_to_virtual_source(const Patient_Parameters_t& pat,
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
    spot.x = pos.x / ((SAD.x + pos.z) / SAD.x);
    spot.y = pos.y / ((SAD.y + pos.z) / SAD.y);
    return spot;
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

