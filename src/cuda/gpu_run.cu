#include "gpu_run.cuh"

#include "gpu_ray_kernel.cuh"
#include "gpu_ray_positioning_kernel.cuh"
#include "gpu_errorcheck.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_utils.cuh"
#include "special_types.hpp"

#include <iostream>
#include <fstream>
#include <numeric>
#include <vector>
#include <sys/stat.h>

void calculateRays(const std::vector<short> spots_per_field,
                   float4* positions_scorer,
                   float4* directions_scorer,
                   short2* metadata_scorer,
                   float* traces_scorer)
{
    short* spots_per_field_gpu;
    array_to_device<short>(spots_per_field_gpu, spots_per_field.data(), spots_per_field.size());
    size_t total_spots = std::accumulate(spots_per_field.begin(), spots_per_field.end(), 0);
    std::cout << std::endl;
    std::cout << "Calculating " << total_spots << " rays ..." << std::endl;
    int nblocks = 1 + (total_spots-1)/NTHREAD_PER_BLOCK_RAYS;
    calculateRays_kernel<<<nblocks, NTHREAD_PER_BLOCK_RAYS>>>(total_spots,
                                                              spots_per_field_gpu,
                                                              positions_scorer,
                                                              directions_scorer,
                                                              metadata_scorer,
                                                              traces_scorer);
    check_kernel_execution(__FILE__, __LINE__);
    cudaFree(spots_per_field_gpu);
}

void buffers_to_device(const std::vector<float4>& xbuffer,
                       const std::vector<float4>& vxbuffer,
                       const std::vector<short2>& ixbuffer)
{
    unsigned int num = xbuffer.size();

    // prepare GPU
    size_t bytes1 = sizeof(float4)*num;
    size_t bytes2 = sizeof(short2)*num;
    gpuErrchk( cudaMalloc((void **) &xdata,  bytes1) );
    gpuErrchk( cudaMalloc((void **) &vxdata, bytes1) );
    gpuErrchk( cudaMalloc((void **) &ixdata, bytes2) );
    gpuErrchk( cudaMemcpyToSymbol(xdata,  xbuffer.data(),  bytes1, 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(vxdata, vxbuffer.data(), bytes1, 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(ixdata, ixbuffer.data(), bytes2, 0, cudaMemcpyHostToDevice) );
}

void virtual_src_to_treatment_plane(const unsigned int num,
                                    const std::vector<BeamAngles_t>& angles,
                                    const float3& ct_offsets)
{
    std::vector<float2> temp(angles.size());
    for (size_t i = 0; i < angles.size(); i++)
    {
        temp[i].x = angles.at(i).gantry;
        temp[i].y = angles.at(i).couch;
    }

    float2* angles_gpu;
    array_to_device<float2>(angles_gpu, temp.data(), angles.size());

    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_SOURCE;
    virtual_src_to_treatment_plane_kernel<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(num, angles_gpu, ct_offsets);
    check_kernel_execution(__FILE__, __LINE__);

    cudaFree(angles_gpu);
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
