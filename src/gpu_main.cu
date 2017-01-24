#include "gpu_main.cuh"

#include "patient_parameters.hpp"
#include "gpu_errorcheck.cuh"
#include "gpu_run.cuh"
#include "gpu_geometry_operations.cuh"
#include "gpu_main.cuh"

#include "tramp.hpp"
#include "gpu_run.cuh"
#include "gpu_device_interaction.cuh"
#include "gpu_ct_to_device.cuh"
#include "volume.hpp"

#include <iostream>
#include <string>
#include <vector>


void gpu_launch(const Patient_Parameters_t &pat, const Patient_Volume_t &ct)
{
    int device = 0;
    cudaSetDevice(device);
    printDevProp(device, false);

    runCalculation(pat, ct);
}

void runCalculation(const Patient_Parameters_t &pat, const Patient_Volume_t &ct)
{
    // mark the start total time timer
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    gpu_ct_to_device::setDimensions(ct);
    gpu_ct_to_device::setDensities(ct);

    // the simulation is initialized once, but the calculation is launched nbeams_h times
    for(size_t i=0; i < pat.nbeams; i++)
    {
        // Create tramp object
        Tramp_t tramp(pat.tramp_files.at(i));
        // Create scorer array
        gpuErrchk( cudaMemcpyToSymbol(nspots, &tramp.nspots, sizeof(unsigned int), 0, cudaMemcpyHostToDevice) );
        gpuErrchk( cudaMalloc( (void **) &scorer, sizeof(float3)*tramp.nspots) );

        std::vector<float4> xbuffer;
        std::vector<float4> vxbuffer;

        init_rays(pat, i, xbuffer, vxbuffer);
        float3 ct_offsets = make_float3(pat.ct.offset.x, pat.ct.offset.y, pat.ct.offset.z);
        calculateRays(xbuffer, vxbuffer, pat.angles.at(i), ct_offsets);
        // outputScorerResults(i);
        clearScorer(scorer, sizeof(float3)*tramp.nspots);
        std::cout << std::endl;
    }

    // Finalize the entire computation
    freeMemory();

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float dt_ms;
    cudaEventElapsedTime(&dt_ms, start, stop);
    cudaThreadExit();
    cudaDeviceReset();

    std::cout << std::endl;
    std::cout << "Program time:  "<< dt_ms << "  (ms)" << std::endl;
    std::cout << "Have a nice day!" << std::endl;
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
