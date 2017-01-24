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
#include <vector>

void init_rays(const Patient_Parameters_t& pat,
               const unsigned int ibeam,
               std::vector<float4>& xbuffer,
               std::vector<float4>& vxbuffer)
//  initialize particle buffer
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
        std::cout << "    Range shifter z upstream:   " << rs.zup   << " cm" << std::endl;
        std::cout << "    Range shifter z downstream: " << rs.zdown << " cm" << std::endl;
        src.zeff = abs(rs.zup) > abs(src.zeff) ? rs.zup : src.zeff;
    }
    if(pat.apertures[ibeam].exists)
    {
        Aperture_Dims_t ap = pat.apertures[ibeam];
        std::cout << "    Aperture thickness:    " << ap.thick << " cm" << std::endl;
        std::cout << "    Aperture z downstream: " << ap.zdown << " cm" << std::endl;
    }
    std::cout << "    Source Z plane: " << src.z << " cm" << std::endl;
    std::cout << "    Effective source Z plane: " << src.zeff << " cm" << std::endl;

    float2 SAD = make_float2(pat.virtualSAD.a, pat.virtualSAD.b); // cm

    xbuffer.reserve(src.nspots);
    vxbuffer.reserve(src.nspots);
    for(unsigned int i=0; i < src.nspots; i++)
    { // LOOP OVER SPOTS
        Spot_t& spot = src.spots[i];
        // From fluence map coordinates to MCAuto coordinates.
        // This is the coordinate system in which the phase space files are written.
        // Then I just have to perform the same transformations as with the phase spaces.
        float3 pos; // cm
        // The x,y coordinates are push together or pulled appart depending on the z coord they are going to be initialized in.
        // Some angular spread will be added too
        pos.x = ((SAD.x + src.zeff) / SAD.x) * spot.x;
        pos.y = ((SAD.y + src.zeff) / SAD.y) * spot.y;
        pos.z = src.zeff;

        float3 dCos;
        float a = (spot.x-pos.x)/abs(src.zeff);
        float b = (spot.y-pos.y)/abs(src.zeff);
        float norm = sqrt(a*a + b*b + 1.f);
        dCos.x = a/norm;
        dCos.y = b/norm;

        float temp = 1.0f - dCos.x*dCos.x - dCos.y*dCos.y;
        if(temp < 0)
        {
            std::cerr << "Something went wrong calculating direction cosines:" << std::endl;
            std::cerr << "Are these correct and they make sense?" << std::endl;
            std::cerr << "    Nominal Z plane:   " << src.z       << std::endl;
            std::cerr << "    Effective Z plane: " << src.zeff    << std::endl;
            std::cerr << "    SAD.x:             " << SAD.x       << std::endl;
            std::cerr << "    SAD.y:             " << SAD.y       << std::endl;
            std::cerr << "    spot #:            " << i           << std::endl;
            std::cerr << "    spot x:            " << spot.x      << std::endl;
            std::cerr << "    spot y:            " << spot.y      << std::endl;
            exit(EXIT_FAILURE);
        };
        dCos.z = sqrt(temp);

        // From energy to WEPL
        float energy = src.energies_internal.at(i)*1000000; // eV
        float wepl   = src.wepls.at(i);

        xbuffer.push_back( make_float4( -pos.y, -pos.x, pos.z, wepl) );
        vxbuffer.push_back( make_float4( -dCos.y, -dCos.x, dCos.z, energy) );
    }
}


void calculateRays(std::vector<float4>& xbuffer,
                   std::vector<float4>& vxbuffer,
                   const BeamAngles_t& ang,
                   const float3& ct_offsets)
{
    float2 angles = make_float2(ang.gantry, ang.couch);

    unsigned int num = rays_to_device(xbuffer, vxbuffer, angles, ct_offsets);

    //      simulate a batch of rays
    std::cout << "\tCalculating " << num << " rays ..." << std::endl;
    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_RAYS;
    calculateRays_kernel<<<nblocks, NTHREAD_PER_BLOCK_RAYS>>>(num, scorer);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaThreadSynchronize() );
}


unsigned int rays_to_device(std::vector<float4>& xbuffer,
                            std::vector<float4>& vxbuffer,
                            const float2& angles,
                            const float3& ct_offsets)
{
    unsigned int total_rays = xbuffer.size();

    unsigned int num = min(NRAYS, total_rays);

    // prepare GPU
    size_t bytes = sizeof(float4)*num;
    gpuErrchk( cudaMalloc((void **) &xdata, bytes) );
    gpuErrchk( cudaMalloc((void **) &vxdata, bytes) );
    gpuErrchk( cudaMemcpyToSymbol(xdata, xbuffer.data(), bytes, 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(vxdata, vxbuffer.data(), bytes, 0, cudaMemcpyHostToDevice) );

    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_SOURCE;
    rays_to_device_kernel<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(num, angles, ct_offsets);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaThreadSynchronize() );

    return num;
}


void clearScorer(void *s, size_t sz)
{
    cudaMemset(s, 0, sz);
}

void clearScorer()
{
    cudaMemset(scorer, 0, sizeof(float3)*nspots);
}

void freeMemory()
{ 
    cudaFreeArray(dens);
    cudaUnbindTexture(dens_tex);
    cudaFree(scorer);
}


