#include "gpu_main.cuh"

#include "gpu_ct_to_device.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_errorcheck.cuh"
#include "gpu_physics_data_to_device.cuh"
#include "gpu_run.cuh"
#include "gpu_source_positioning.cuh"
#include "gpu_utils.cuh"
#include "initialize_rays.cuh"
#include "patient_parameters.hpp"
#include "utils.hpp"
#include "volume.hpp"

#include <iostream>
#include <string>
#include <vector>

void gpu_raytrace_original (const Patient_Parameters_t& pat,
                            const Volume_t& ct,
                            Array4<float>& endpoints,
                            Array4<float>& init_pos,
                            Array4<float>& init_pat_pos,
                            std::string output_file)
{
    // Set geometry in GPU
    gpu_ct_to_device::sendGeometries(ct);

    // Create host buffers and initialize rays
    std::vector<float4> xbuffer;
    std::vector<float4> vxbuffer;
    std::vector<short2> ixbuffer;
    create_virtual_source_buffers (pat, xbuffer, vxbuffer, ixbuffer);
    buffers_to_device (xbuffer, vxbuffer, ixbuffer);
    virtual_src_to_treatment_plane (xbuffer.size(), pat.angles,
                                    make_float3(pat.ct.offset.x, pat.ct.offset.y, pat.ct.offset.z));

    gpuErrchk( cudaMemcpyFromSymbol(&(init_pat_pos[0].x), xdata, sizeof(float4)*xbuffer.size(), 0, cudaMemcpyDeviceToHost) );
    // Copy buffer with initial positions and wepl
    for (size_t i = 0; i < xbuffer.size(); i++)
    {
        init_pos.at(i).x = xbuffer[i].x;
        init_pos.at(i).y = xbuffer[i].y;
        init_pos.at(i).z = xbuffer[i].z;
        init_pos.at(i).w = xbuffer[i].w; // wepl
    }

    gpu_raytrace (pat, endpoints, output_file);
}

void gpu_raytrace_warped (const Patient_Parameters_t &pat,
                          const Volume_t &ct,
                          const Array4<float>& orig_endpoints,
                          const Array4<float>& init_pos,
                          Array4<float>& endpoints,
                          std::string output_file)
{
    // Set geometry in GPU
    gpu_ct_to_device::sendGeometries(ct);

    // Create host buffers and initialize rays
    std::vector<float4> xbuffer;
    std::vector<float4> vxbuffer;
    std::vector<short2> ixbuffer;
    create_treatment_plane_buffers (pat, orig_endpoints, init_pos,
                                    xbuffer, vxbuffer, ixbuffer);
    buffers_to_device (xbuffer, vxbuffer, ixbuffer);

    gpu_raytrace (pat, endpoints, output_file, orig_endpoints);
}

void gpu_raytrace (const Patient_Parameters_t& pat,
                   Array4<float>& endpoints,
                   std::string output_file,
                   const Array4<float>& orig_endpoints)
{
    // Create scorer array
    float4* pos_scorer = NULL;
    allocate_scorer<float4>(pos_scorer, pat.total_spots);

    // Calculate rays
    if (output_file.empty())
    {
        do_raytrace (pat.spots_per_field, pos_scorer, NULL, orig_endpoints);
    }
    else
    {
        float* traces_scorer = NULL;
        allocate_scorer<float>(traces_scorer, pat.ct.total);
        do_raytrace (pat.spots_per_field, pos_scorer, traces_scorer, orig_endpoints);
        Volume_t traces(pat.ct);
        retrieve_scorer<float, float>(&traces.data[0], traces_scorer, traces.nElements);
        // traces.output("output_volume.mha", "mha");
        traces.output(output_file, "bin");
        gpuErrchk( cudaFree(traces_scorer) );
    }

    retrieve_scorer<float, float4>(&(endpoints[0].x), pos_scorer, pat.total_spots);
    // Free memory
    gpuErrchk( cudaFree(pos_scorer) );
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

    gpu_physics_to_device::sendWaterRestrictedSPower();
    gpu_physics_to_device::sendMassStoppingPowerRatio();
}

void stop_device(cudaEvent_t& start)
{
    freePhysicsMemory();

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
