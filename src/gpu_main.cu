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

void initialize_device(cudaEvent_t& start, cudaEvent_t& stop)
{
	// mark the start total time timer
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start);

	// Set device
	int device = 0;
	cudaSetDevice(device);
	bool verbose = true;
	printDevProp(device, verbose);
}

void stop_device(cudaEvent_t& start, cudaEvent_t& stop)
{
	// Get timing
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float dt_ms;
	cudaEventElapsedTime(&dt_ms, start, stop);
	cudaThreadExit();
	cudaDeviceReset();

	std::cout << std::endl;
	std::cout << "Tracing time: "  << dt_ms/1000 << " s" << std::endl;
}

std::vector<float4> gpu_get_beam_endpoints(const Patient_Parameters_t &pat,
                                           const Patient_Volume_t &ct)
{
    // Run
    std::vector<float4> endpoints(pat.total_spots);
    runCalculation(pat, ct, endpoints);
    


    return endpoints;
}

void runCalculation(const Patient_Parameters_t &pat,
                    const Patient_Volume_t &ct,
                    std::vector<float4>& endpoints)
{
    // Set geometry in GPU
    gpu_ct_to_device::setDimensions(ct);
    gpu_ct_to_device::setDensities(ct); // density correction file
                                        // is an additional parameter
                                        // with a default value

    // Create scorer array
    float4* endpoints_scorer = NULL;
    gpuErrchk( cudaMalloc( (void **) &endpoints_scorer, sizeof(float4)*pat.total_spots) );
    gpuErrchk( cudaMemset( (void *) endpoints_scorer, 0, sizeof(float4)*pat.total_spots) );
    float* traces_scorer = NULL;
    gpuErrchk( cudaMalloc( (void **) &traces_scorer, sizeof(float)*pat.ct.total) );
    gpuErrchk( cudaMemset( (void *) traces_scorer, 0, sizeof(float)*pat.ct.total) );

    // Create host buffers and fill them
    std::vector<float4> xbuffer;
    std::vector<float4> vxbuffer;
    std::vector<short2> ixbuffer;
    init_rays(pat, xbuffer, vxbuffer, ixbuffer);

    // Get ct offsets
    float3 ct_offsets = make_float3(pat.ct.offset.x, pat.ct.offset.y, pat.ct.offset.z);

    // Calculate rays
    calculateRays(xbuffer, vxbuffer, ixbuffer,
                  pat.angles, pat.spots_per_field.data(),
                  ct_offsets, endpoints_scorer, traces_scorer);

    gpuErrchk( cudaMemcpy(&endpoints[0], endpoints_scorer, sizeof(float4)*pat.total_spots, cudaMemcpyDeviceToHost) );

    Patient_Volume_t traces(pat.ct);
    gpuErrchk( cudaMemcpy(&traces.hu[0], traces_scorer, sizeof(float)*traces.nElements, cudaMemcpyDeviceToHost) );
    traces.output("output_volume.bin");

    // Free memory
    gpuErrchk( cudaFree(endpoints_scorer) );
    gpuErrchk( cudaFree(traces_scorer) );
    freeCTMemory();
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
