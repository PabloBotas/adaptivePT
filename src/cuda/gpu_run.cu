#include "gpu_device_interaction.cuh"
#include "gpu_run.cuh"
#include "gpu_ray_kernel.cuh"
#include "gpu_errorcheck.cuh"
#include "gpu_geometry_operations.cuh"

#include "patient_parameters.hpp"
#include "special_types.hpp"
#include "tramp.hpp"
#include "spot.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>

void init_rays(const Patient_Parameters_t& pat,
               std::vector<float4>& xbuffer,
               std::vector<float4>& vxbuffer,
               std::vector<short2>& ixbuffer)
//  initialize particle buffer
{
    for (unsigned int ibeam = 0; ibeam < pat.beam_names.size(); ibeam++)
    {
        std::cout << "Loading beam: " << pat.beam_names.at(ibeam) << std::endl;
        std::cout << "    Beam number:  " << ibeam << std::endl;
        std::cout << "    Machine:      " << pat.machine << std::endl;
        std::cout << "    Tramp file:   " << pat.tramp_files.at(ibeam) << std::endl;

        Tramp_t src(pat.tramp_files.at(ibeam), pat.machine);
        src.z = pat.isocenter_to_beam_distance.at(ibeam); // cm
        src.zeff = src.z;

        if(pat.range_shifters[ibeam].exists)
        {
            RangeShifter_Dims_t rs = pat.range_shifters[ibeam];
            std::cout << "    Range shifter thickness:    " << rs.thick << " cm" << std::endl;
            for (size_t i = 0; i < src.wepls.size(); i++)
            {
                src.wepls.at(i) -= rs.wepl;
            }
        }
        if(pat.apertures[ibeam].exists)
        {
            Aperture_Dims_t ap = pat.apertures[ibeam];
            std::cout << "    Aperture thickness:    " << ap.thick << " cm" << std::endl;
            std::cout << "    Aperture z downstream: " << ap.zdown << " cm" << std::endl;
        }
        std::cout << "    Source Z plane: " << src.z << " cm" << std::endl;

        float2 SAD = make_float2(pat.virtualSAD.a, pat.virtualSAD.b); // cm

        xbuffer.reserve(src.nspots);
        vxbuffer.reserve(src.nspots);
        ixbuffer.reserve(src.nspots);
        for(unsigned int i=0; i < src.nspots; i++)
        { // LOOP OVER SPOTS
            Spot_t& spot = src.spots[i];

            float3 pos; // cm
            pos.x = ((SAD.x + src.z) / SAD.x) * spot.x;
            pos.y = ((SAD.y + src.z) / SAD.y) * spot.y;
            pos.z = src.z;

            float3 dCos;
            float a = (spot.x-pos.x)/abs(src.z);
            float b = (spot.y-pos.y)/abs(src.z);
            float norm = sqrt(a*a + b*b + 1.f);
            dCos.x = a/norm;
            dCos.y = b/norm;

            float temp = 1.0f - dCos.x*dCos.x - dCos.y*dCos.y;
            if(temp < 0)
            {
                std::cerr << "Something went wrong calculating direction cosines:" << std::endl;
                std::cerr << "    Z plane: " << src.z  << std::endl;
                std::cerr << "    SAD.x:   " << SAD.x  << std::endl;
                std::cerr << "    SAD.y:   " << SAD.y  << std::endl;
                std::cerr << "    spot #:  " << i      << std::endl;
                std::cerr << "    spot x:  " << spot.x << std::endl;
                std::cerr << "    spot y:  " << spot.y << std::endl;
                exit(EXIT_FAILURE);
            };
            dCos.z = sqrt(temp);

            // From energy to WEPL
            float energy = src.energies_internal.at(i)*1000000; // eV
            float wepl   = src.wepls.at(i);                     // cm

            xbuffer.push_back( make_float4(-pos.y, -pos.x, pos.z, wepl) );
            vxbuffer.push_back( make_float4(-dCos.y, -dCos.x, dCos.z, energy) );
            ixbuffer.push_back( make_short2(ibeam, i) );
        }
    }
}


void calculateRays(const std::vector<float4>& xbuffer,
                   const std::vector<float4>& vxbuffer,
                   const std::vector<short2>& ixbuffer,
                   const std::vector<BeamAngles_t>& angles,
                   const short* spots_per_beam,
                   const float3& ct_offsets,
                   float4* endpoints_scorer,
                   float* traces_scorer)
{
    unsigned int total_spots = rays_to_device(xbuffer, vxbuffer, ixbuffer, angles, ct_offsets);

    short* spots_per_beam_gpu;
    gpuErrchk( cudaMalloc((void **) &spots_per_beam_gpu, sizeof(short)*angles.size()) );
    gpuErrchk( cudaMemcpy(spots_per_beam_gpu, spots_per_beam, sizeof(short)*angles.size(), cudaMemcpyHostToDevice) );
    //      simulate a batch of rays
    std::cout << std::endl;
    std::cout << "Calculating " << total_spots << " rays ..." << std::endl;
    int nblocks = 1 + (total_spots-1)/NTHREAD_PER_BLOCK_RAYS;
    calculateRays_kernel<<<nblocks, NTHREAD_PER_BLOCK_RAYS>>>(total_spots, endpoints_scorer, traces_scorer, spots_per_beam_gpu);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaThreadSynchronize() );
}

unsigned int rays_to_device(const std::vector<float4>& xbuffer,
                            const std::vector<float4>& vxbuffer,
                            const std::vector<short2>& ixbuffer,
                            const std::vector<BeamAngles_t>& angles,
                            const float3& ct_offsets)
{
    unsigned int num = xbuffer.size();

    // prepare GPU
    size_t bytes1 = sizeof(float4)*num;
    size_t bytes2 = sizeof(short2)*num;
    gpuErrchk( cudaMalloc((void **) &xdata, bytes1) );
    gpuErrchk( cudaMalloc((void **) &vxdata, bytes1) );
    gpuErrchk( cudaMalloc((void **) &ixdata, bytes2) );
    gpuErrchk( cudaMemcpyToSymbol(xdata,  xbuffer.data(), bytes1, 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(vxdata, vxbuffer.data(), bytes1, 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(ixdata, ixbuffer.data(), bytes2, 0, cudaMemcpyHostToDevice) );

    std::vector<float2> temp(angles.size());
    for (size_t i = 0; i < angles.size(); i++)
    {
        temp[i].x = angles.at(i).gantry;
        temp[i].y = angles.at(i).couch;
    }

    float2* angles_gpu;
    gpuErrchk( cudaMalloc((void **) &angles_gpu, sizeof(float2)*angles.size()) );
    gpuErrchk( cudaMemcpy(angles_gpu, temp.data(), sizeof(float2)*angles.size(), cudaMemcpyHostToDevice) );

    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_SOURCE;
    rays_to_device_kernel<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(num, angles_gpu, ct_offsets);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaThreadSynchronize() );

    cudaFree(angles_gpu);

    return num;
}

void freeCTMemory()
{
    cudaFreeArray(dens);
    cudaUnbindTexture(dens_tex);
}
