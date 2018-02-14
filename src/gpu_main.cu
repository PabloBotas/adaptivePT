#include "gpu_main.cuh"

#include "command_line_parser.hpp"
#include "gpu_ct_to_device.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_errorcheck.cuh"
#include "gpu_influence_kernel.cuh"
#include "gpu_physics_data_to_device.cuh"
#include "gpu_run.cuh"
#include "gpu_source_positioning.cuh"
#include "gpu_utils.cuh"
#include "initialize_rays.cuh"
#include "patient_parameters.hpp"
#include "tramp.hpp"
#include "utils.hpp"
#include "volume.hpp"

#include <iostream>
#include <string>
#include <vector>

void gpu_raytrace_original (const Patient_Parameters_t& pat,
                            const Parser& parser,
                            const Volume_t& ct,
                            Array4<float>& endpoints,
                            Array4<float>& initpos_xbuffer_dbg,
                            Array4<float>& initpos)
{
    // Set geometry in GPU
    gpu_ct_to_device::sendGeometries(ct);
    gpu_ct_to_device::sendMask(parser.v_field_mask_files, pat.ct);

    // Create host buffers and initialize rays
    std::vector<float4> xbuffer, vxbuffer;
    std::vector<short2> ixbuffer;
    create_virtual_source_buffers (pat, xbuffer, vxbuffer, ixbuffer);
    buffers_to_device (xbuffer, vxbuffer, ixbuffer, true);
    virtual_src_to_treatment_plane (xbuffer.size(), pat.angles,
                                    make_float3(pat.ct.offset.x, pat.ct.offset.y, pat.ct.offset.z));

    gpuErrchk( cudaMemcpyFromSymbol(&(initpos[0].x), xdata, sizeof(float4)*xbuffer.size(), 0, cudaMemcpyDeviceToHost) );
    // Copy buffer with initial positions and energy
    for (size_t i = 0; i < xbuffer.size(); i++) {
        initpos_xbuffer_dbg.at(i).x = xbuffer[i].x;
        initpos_xbuffer_dbg.at(i).y = xbuffer[i].y;
        initpos_xbuffer_dbg.at(i).z = xbuffer[i].z;
        initpos_xbuffer_dbg.at(i).w = xbuffer[i].w; // energy
    }

    Array3<float> dummy_not_used;
    gpu_raytrace (pat, ct, endpoints, parser.ct_traces_file, dummy_not_used);

    gpu_ct_to_device::removeMask();
}

void gpu_raytrace_warped (const Patient_Parameters_t &pat,     // pat.ct is CBCT, original_ct is CT
                          const Parser& parser,
                          const Volume_t &cbct,                // CBCT
                          const Array4<float>& orig_endpoints, // input positions in CT space
                          const Array4<float>& init_pos,       // input positions in CT space
                          Array4<float>& endpoints,
                          std::vector<float>& new_energies,
                          Array3<float>& traces_info)
{
    // Set geometry in GPU
    gpu_ct_to_device::sendGeometries(cbct);

    // Create host buffers and initialize rays
    std::vector<float4> xbuffer;
    std::vector<float4> vxbuffer;
    std::vector<short2> ixbuffer;
    set_treatment_plane_buffers_ct_space (pat, orig_endpoints, init_pos,
                                          xbuffer, vxbuffer, ixbuffer);
    bool need_reallocation = false;
    buffers_to_device (xbuffer, vxbuffer, ixbuffer, need_reallocation);

    correct_offsets (xbuffer.size(),
                    make_float3(cbct.origin.x, cbct.origin.y, cbct.origin.z),
                    make_int3(cbct.n.x, cbct.n.y, cbct.n.z),
                    make_float3(cbct.d.x, cbct.d.y, cbct.d.z),
                    make_float3(pat.ct.origin.x, pat.ct.origin.y, pat.ct.origin.z),
                    make_int3(pat.ct.n.x, pat.ct.n.y, pat.ct.n.z),
                    make_float3(pat.ct.d.x, pat.ct.d.y, pat.ct.d.z));
    Array4<float> off_endpoints = offset_endpoints (orig_endpoints, 
                    make_float3(cbct.origin.x, cbct.origin.y, cbct.origin.z),
                    make_int3(cbct.n.x, cbct.n.y, cbct.n.z),
                    make_float3(cbct.d.x, cbct.d.y, cbct.d.z),
                    make_float3(pat.ct.origin.x, pat.ct.origin.y, pat.ct.origin.z),
                    make_int3(pat.ct.n.x, pat.ct.n.y, pat.ct.n.z),
                    make_float3(pat.ct.d.x, pat.ct.d.y, pat.ct.d.z));

    gpu_raytrace (pat, cbct, endpoints, parser.cbct_traces_file, traces_info, off_endpoints);

    new_energies.resize(endpoints.size());
    for (size_t i = 0; i < new_energies.size(); ++i) {
        new_energies.at(i) = endpoints.at(i).w;
    }
}

void gpu_raytrace (const Patient_Parameters_t& pat,
                   const Volume_t &vol,
                   Array4<float>& endpoints,
                   std::string traces_file,
                   Array3<float>& traces_info,
                   const Array4<float>& orig_endpoints)
{
    // Create scorer array
    float4* pos_scorer = NULL;
    allocate_scorer<float4>(pos_scorer, pat.total_spots);

    // Calculate rays
    if (traces_file.empty()) {
        do_raytrace(pat.spots_per_field, pos_scorer, NULL, orig_endpoints, traces_info);
    } else {
        float* traces_scorer = NULL;
        allocate_scorer<float>(traces_scorer, vol.nElements);
        do_raytrace(pat.spots_per_field, pos_scorer, traces_scorer, orig_endpoints, traces_info);
        Volume_t traces(vol.getMetadata());
        retrieve_scorer<float, float>(&traces.data[0], traces_scorer, traces.nElements);
        // traces.output("output_volume.mha");
        traces.output(traces_file);
        gpuErrchk( cudaFree(traces_scorer) );
    }

    retrieve_scorer<float, float4>(&(endpoints[0].x), pos_scorer, pat.total_spots);
    // Free memory
    gpuErrchk( cudaFree(pos_scorer) );
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
    // cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 1024000*30);

    gpu_physics_to_device::sendWaterRestrictedSPower();
    gpu_physics_to_device::sendMassStoppingPowerRatio();
    gpu_physics_to_device::sendBraggPeakFits();
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
    if (verbose) {
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
    } else {
        std::cout << "Using card (" << device << "): " << devProp.name << std::endl;
    }
}
