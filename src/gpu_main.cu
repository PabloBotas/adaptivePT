#include "gpu_main.cuh"

#include "patient_parameters.hpp"
#include "volume.hpp"
#include "gpu_ct_to_device.cuh"
#include "initialize_rays.cuh"
#include "gpu_errorcheck.cuh"
#include "gpu_run.cuh"
#include "gpu_source_positioning.cuh"
#include "gpu_utils.cuh"
#include "utils.hpp"

#include <iostream>
#include <string>
#include <vector>

std::vector< Vector4_t<float> > gpu_raytrace_plan(const Patient_Parameters_t &pat,
                                                  const Patient_Volume_t &ct)
{
    // Run
    std::vector< Vector4_t<float> > endpoints(pat.total_spots);
    gpu_raytrace_plan(pat, ct, endpoints);
    // utils::cm_to_mm(endpoints);

    return endpoints;
}

void gpu_raytrace_plan(const Patient_Parameters_t& pat,
                       const Patient_Volume_t& ct,
                       std::vector< Vector4_t<float> >& endpoints)
{
    // Set geometry in GPU
    gpu_ct_to_device::sendDimensions(ct);
    std::vector<int> HU_indexes = gpu_ct_to_device::sendMassStoppingPowerRatio();
    gpu_ct_to_device::sendDensities(ct);
    gpu_ct_to_device::sendMaterialId(ct, HU_indexes);

    // Create host buffers and initialize rays
    std::vector<float4> xbuffer;
    std::vector<float4> vxbuffer;
    std::vector<short2> ixbuffer;
    create_virtual_source_buffers(pat, xbuffer, vxbuffer, ixbuffer);
    buffers_to_device(xbuffer, vxbuffer, ixbuffer);
    virtual_src_to_treatment_plane(xbuffer.size(), pat.angles,
                                   make_float3(pat.ct.offset.x, pat.ct.offset.y, pat.ct.offset.z));

    // Create scorer array
    float4* pos_scorer = NULL;
    float4* dir_scorer = NULL;
    short2* meta_scorer = NULL;
    allocate_scorer<float4>(pos_scorer, pat.total_spots);
    allocate_scorer<float4>(dir_scorer, pat.total_spots);
    allocate_scorer<short2>(meta_scorer, pat.total_spots);

    // Calculate rays
#if !defined __OUTPUT_SCORER_VOLUME__
    calculateRays(pat.spots_per_field,
                  pos_scorer, dir_scorer, meta_scorer);
#else
    float* traces_scorer = NULL;
    allocate_scorer<float>(traces_scorer, pat.ct.total);
    calculateRays(pat.spots_per_field,
                  pos_scorer, dir_scorer, meta_scorer,
                  traces_scorer);
    Patient_Volume_t traces(pat.ct);
    gpuErrchk( cudaMemcpy(&traces.hu[0], traces_scorer, sizeof(float)*traces.nElements, cudaMemcpyDeviceToHost) );
    // traces.output("output_volume.mha", "mha");
    traces.output("output_volume.raw", "bin");
    gpuErrchk( cudaFree(traces_scorer) );
#endif
    retrieve_scorer<float, float4>(&(endpoints[0].x), pos_scorer, pat.total_spots);

    // Free memory
    gpuErrchk( cudaFree(pos_scorer) );
    gpuErrchk( cudaFree(dir_scorer) );
    gpuErrchk( cudaFree(meta_scorer) );
    freeCTMemory();
}

void initialize_device(cudaEvent_t& start)
{
    // mark the start total time timer
    cudaEventCreate(&start);
    cudaEventRecord(start);

    // Set device
    int device = 0;
    cudaSetDevice(device);
    bool verbose = false;
    printDevProp(device, verbose);
}

void stop_device(cudaEvent_t& start)
{
    // Get timing
    cudaEvent_t stop;
    cudaEventCreate(&stop);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float dt_ms;
    cudaEventElapsedTime(&dt_ms, start, stop);
    cudaThreadExit();
    cudaDeviceReset();

    std::cout << std::endl;
    std::cout << "Tracing time: "  << dt_ms/1000 << " s" << std::endl;
}

void printDevProp(const int device, bool verbose)
{
    cudaDeviceProp devProp;
    cudaGetDeviceProperties(&devProp, device);
    if(verbose)
    {
        std::cout << "Using device #:        " << device << std::endl;
        std::cout << "Name:                  " << devProp.name << std::endl;
        std::cout << "Compute capability:    " << devProp.major << "." << devProp.minor << std::endl;
        std::cout << "Global memory:         " << devProp.totalGlobalMem/1024.0/1024.0 << " MB" << std::endl;
        std::cout << "Shared memory /block:  " << devProp.sharedMemPerBlock/1024.0 << std::endl;
        std::cout << "Registers /block:      " << devProp.regsPerBlock << std::endl;
        std::cout << "Warp size:             " << devProp.warpSize << std::endl;
        std::cout << "Memory pitch:          " << devProp.memPitch << std::endl;
        std::cout << "Threads /block:        " << devProp.maxThreadsPerBlock << std::endl;
        std::cout << "Maximum dim of block:  " << devProp.maxThreadsDim[0] << "," << devProp.maxThreadsDim[1] << "," << devProp.maxThreadsDim[2] << std::endl;
        std::cout << "Maximum dim of grid:   " << devProp.maxGridSize[0] << "," << devProp.maxGridSize[1] << "," << devProp.maxGridSize[2] << std::endl;
        std::cout << "Clock rate:            " << devProp.clockRate/1000000.0 << " GHz" << std::endl;
        std::cout << "Total constant memory: " << devProp.totalConstMem/1024.0 << std::endl;
        std::cout << "Texture alignment:     " << devProp.textureAlignment << std::endl;
        std::cout << "Concurrent copy/exec:  " << (devProp.deviceOverlap ? "Yes" : "No") << std::endl;
        std::cout << "Multiprocessors:       " << devProp.multiProcessorCount << std::endl;
        std::cout << "Kernel timeout:        " << (devProp.kernelExecTimeoutEnabled ? "Yes" : "No") << std::endl;
    }
    else
    {
        std::cout << "Using card (" << device << "): " << devProp.name << std::endl;
    }
}
