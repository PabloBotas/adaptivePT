#include "gpu_run.cuh"

#include "gpu_ray_kernel.cuh"
#include "gpu_ray_positioning_kernel.cuh"
#include "gpu_errorcheck.cuh"
#include "gpu_device_globals.cuh"
#include "special_types.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>

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
    calculateRays_kernel<<<nblocks, NTHREAD_PER_BLOCK_RAYS>>>(total_spots, spots_per_beam_gpu, endpoints_scorer, traces_scorer);
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
    rays_to_delivery_plane<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(num, angles_gpu, ct_offsets);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaThreadSynchronize() );

    cudaFree(angles_gpu);

    return num;
}

void freeCTMemory()
{
    cudaFreeArray(dens);
    cudaUnbindTexture(dens_tex);
    cudaFreeArray(matid);
    cudaUnbindTexture(matid_tex);
    cudaFreeArray(stp_ratio_array);
    cudaUnbindTexture(stp_ratio_tex);
}
