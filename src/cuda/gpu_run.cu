#include "gpu_run.cuh"

#include "gpu_ray_kernel.cuh"
#include "gpu_source_positioning.cuh"
#include "gpu_errorcheck.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_utils.cuh"
#include "special_types.hpp"

#include <iostream>
#include <fstream>
#include <numeric>
#include <vector>
#include <sys/stat.h>

void do_raytrace (const std::vector<short>& spots_per_field,
                  float4* positions_scorer,
                  float* traces_scorer,
                  const Array4<float>& orig_endpoints,
                  Array3<float>& traces_data)
{
    ushort total_spots = std::accumulate(spots_per_field.begin(), spots_per_field.end(), 0);
    // Set up optional target endpoints
    float4* dev_orig_endpoints = NULL;
    float3* dev_traces_data = NULL;
    if (!orig_endpoints.empty()) {
        array_to_device<float4, Vector4_t<float>>(dev_orig_endpoints,
                                                  orig_endpoints.data(), orig_endpoints.size());
        // Allocate for start-end position per spot
        allocate_scorer<float3>(dev_traces_data, total_spots*2);
    }

    short* spf_gpu = NULL;
    array_to_device<short>(spf_gpu, spots_per_field.data(), spots_per_field.size());
    
    std::cout << std::endl;
    std::cout << "Calculating " << total_spots << " rays ..." << std::endl;
    int nblocks = 1 + (total_spots-1)/NTHREAD_PER_BLOCK_RAYS;
    raytrace_plan_kernel<<<nblocks, NTHREAD_PER_BLOCK_RAYS>>>(total_spots,
                                                              spf_gpu,
                                                              dev_orig_endpoints,
                                                              positions_scorer,
                                                              traces_scorer,
                                                              dev_traces_data);
    check_kernel_execution(__FILE__, __LINE__);

    if (!orig_endpoints.empty()) {
        traces_data.resize(total_spots*2);
        retrieve_scorer<float, float3>(&(traces_data[0].x), dev_traces_data, total_spots*2);
    }

    cudaFree(spf_gpu);
    cudaFree(dev_orig_endpoints);
    cudaFree(dev_traces_data);
}


void buffers_to_device(const std::vector<float4>& xbuffer,
                       const std::vector<float4>& vxbuffer,
                       const std::vector<short2>& ixbuffer,
                       const bool alloc)
{
    unsigned int num = xbuffer.size();

    // prepare GPU
    size_t bytes1 = sizeof(float4)*num;
    size_t bytes2 = sizeof(short2)*num;
    if (alloc) {
        gpuErrchk( cudaMalloc((void **) &xdata,  bytes1) );
        gpuErrchk( cudaMalloc((void **) &vxdata, bytes1) );
        gpuErrchk( cudaMalloc((void **) &ixdata, bytes2) );
    }
    gpuErrchk( cudaMemcpyToSymbol(xdata,  xbuffer.data(),  bytes1, 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(vxdata, vxbuffer.data(), bytes1, 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(ixdata, ixbuffer.data(), bytes2, 0, cudaMemcpyHostToDevice) );
}

