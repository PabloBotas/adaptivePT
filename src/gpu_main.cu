#include "gpu_main.cuh"

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
                            const Volume_t& ct,
                            Array4<double>& endpoints,
                            Array4<double>& initpos_xbuffer_dbg,
                            Array4<double>& initpos,
                            std::string output_file,
                            Array4<double>& influence)
{
    // Set geometry in GPU
    gpu_ct_to_device::sendGeometries(ct);

    // Create host buffers and initialize rays
    std::vector<double4> xbuffer;
    std::vector<double4> vxbuffer;
    std::vector<short2> ixbuffer;
    create_virtual_source_buffers (pat, xbuffer, vxbuffer, ixbuffer);
    buffers_to_device (xbuffer, vxbuffer, ixbuffer, true);
    virtual_src_to_treatment_plane (xbuffer.size(), pat.angles,
                                    make_double3(pat.ct.offset.x, pat.ct.offset.y, pat.ct.offset.z));

    gpuErrchk( cudaMemcpyFromSymbol(&(initpos[0].x), xdata, sizeof(double4)*xbuffer.size(), 0, cudaMemcpyDeviceToHost) );
    // Copy buffer with initial positions and wepl
    for (size_t i = 0; i < xbuffer.size(); i++)
    {
        initpos_xbuffer_dbg.at(i).x = xbuffer[i].x;
        initpos_xbuffer_dbg.at(i).y = xbuffer[i].y;
        initpos_xbuffer_dbg.at(i).z = xbuffer[i].z;
        initpos_xbuffer_dbg.at(i).w = xbuffer[i].w; // wepl
    }

    gpu_raytrace (pat, endpoints, output_file);

    // INFLUENCE
    std::vector<float> inf_cube(ct.nElements);
    std::vector<float> spot_weights;
    for (size_t i = 0; i < pat.tramp_files.size(); ++i)
    {
        Tramp_t tramp(pat.tramp_files.at(i), pat.machine);
        std::vector<float> w = tramp.get_weights();
        spot_weights.reserve(spot_weights.size() + tramp.nspots);
        spot_weights.insert(spot_weights.end(), w.begin(), w.end());
    }
    gpu_calculate_influence (pat.total_spots, endpoints, influence, spot_weights, inf_cube);
#ifdef __INFLUENCE_MATRICES__
    std::cout << "Writting influence_matrix_CT.dat ..." << std::endl;
    std::ofstream fout("influence_matrix_CT.dat", std::ios::out | std::ios::binary);
    for (size_t i = 0; i < influence.size(); ++i)
        fout.write((char*)&influence[i].w, sizeof(double));

    std::cout << "Writting influence_cube_CT.dat ..." << std::endl;
    std::ofstream fout2("influence_cube_CT.dat", std::ios::out | std::ios::binary);
    fout2.write((char*)&inf_cube[0], inf_cube.size()*sizeof(float));
#endif
    // for (size_t i = 0; i < 10; i++)
    // {
    //     std::cerr << influence[i].x << " " << influence[i].y << " " << influence[i].z << " " << influence[i].w << "\n";
    // }
    freeCTMemory();
}

void gpu_raytrace_warped (const Patient_Parameters_t &pat,
                          const Volume_t &ct,
                          const Array4<double>& orig_endpoints,
                          const Array4<double>& init_pos,
                          Array4<double>& endpoints,
                          std::string output_file,
                          Array4<double>& influence)
{
    // Set geometry in GPU
    gpu_ct_to_device::sendGeometries(ct);

    // Create host buffers and initialize rays
    std::vector<double4> xbuffer;
    std::vector<double4> vxbuffer;
    std::vector<short2> ixbuffer;
    create_treatment_plane_buffers (pat, orig_endpoints, init_pos,
                                    xbuffer, vxbuffer, ixbuffer);
    buffers_to_device (xbuffer, vxbuffer, ixbuffer, false);
    correct_offsets (xbuffer.size(), 
        make_double3(pat.ct.offset.x, pat.ct.offset.y, pat.ct.offset.z),
        make_double3(pat.original_ct.offset.x, pat.original_ct.offset.y, pat.original_ct.offset.z));
    Array4<double> off_endpoints = offset_endpoints (orig_endpoints, 
                                       make_double3(pat.ct.offset.x, pat.ct.offset.y, pat.ct.offset.z),
                                       make_double3(pat.original_ct.offset.x, pat.original_ct.offset.y, pat.original_ct.offset.z));

    gpu_raytrace (pat, endpoints, output_file, off_endpoints);

    std::vector<float> inf_cube(ct.nElements);
    std::vector<float> spot_weights;
    for (size_t i = 0; i < pat.tramp_files.size(); ++i)
    {
        Tramp_t tramp(pat.tramp_files.at(i), pat.machine);
        std::vector<float> w = tramp.get_weights();
        spot_weights.reserve(spot_weights.size() + tramp.nspots);
        spot_weights.insert(spot_weights.end(), w.begin(), w.end());
    }
    gpu_calculate_influence (pat.total_spots, endpoints, influence, spot_weights, inf_cube);
#ifdef __INFLUENCE_MATRICES__
    std::cout << "Writting influence_CBCT.dat ..." << std::endl;
    std::ofstream fout("influence_matrix_CBCT.dat", std::ios::out | std::ios::binary);
    for (size_t i = 0; i < influence.size(); ++i)
        fout.write((char*)&influence[i].w, sizeof(double));
    
    std::cout << "Writting influence_cube_CBCT.dat ..." << std::endl;
    std::ofstream fout2("influence_cube_CBCT.dat", std::ios::out | std::ios::binary);
    fout2.write((char*)&inf_cube[0], inf_cube.size()*sizeof(float));
#endif
    // for (size_t i = 0; i < influence.size(); i+=6421)
    // {
    //     std::cout << influence[i].x << " " << influence[i].y << " " << influence[i].z << " " << influence[i].w << "\n";
    // }
    freeCTMemory();
}

void gpu_raytrace (const Patient_Parameters_t& pat,
                   Array4<double>& endpoints,
                   std::string output_file,
                   const Array4<double>& orig_endpoints)
{
    // Create scorer array
    double4* pos_scorer = NULL;
    allocate_scorer<double4>(pos_scorer, pat.total_spots);

    // Calculate rays
    if (output_file.empty())
    {
        do_raytrace(pat.spots_per_field, pos_scorer, NULL, orig_endpoints);
    }
    else
    {
        float* traces_scorer = NULL;
        allocate_scorer<float>(traces_scorer, pat.ct.total);
        do_raytrace(pat.spots_per_field, pos_scorer, traces_scorer, orig_endpoints);
        bool long_data = false;
        Volume_t traces(pat.ct, long_data);
        retrieve_scorer<float, float>(&traces.data[0], traces_scorer, traces.nElements);
        // traces.output("output_volume.mha");
        traces.output(output_file);
        gpuErrchk( cudaFree(traces_scorer) );
    }

    retrieve_scorer<double, double4>(&(endpoints[0].x), pos_scorer, pat.total_spots);
    // Free memory
    gpuErrchk( cudaFree(pos_scorer) );
}

void gpu_calculate_influence (const short& nspots,
                              const Array4<double>& endpoints,
                              Array4<double>& influence,
                              std::vector<float>& spot_weights,
                              std::vector<float>& inf_volume)
{
    // Copy data to device
    double4 *dev_endpoints = NULL;
    array_to_device<double4, Vector4_t<double> >(dev_endpoints, endpoints.data(), endpoints.size());
    // Create scorer array
    double4 *dev_influence = NULL;
    allocate_scorer<double4>(dev_influence, nspots*nspots);
    // Create scorer array
    float *dev_inf_volume = NULL;
    allocate_scorer<float>(dev_inf_volume, inf_volume.size());
    // Create weights array
    float *dev_spot_weights = NULL;
    array_to_device<float>(dev_spot_weights, &spot_weights[0], spot_weights.size());
    
    // Launch influence kernel
    std::cout << "Calculating influence matrix ..." << std::endl;

    int nblocks = 1 + (nspots*nspots-1)/NTHREAD_PER_BLOCK_INFLUENCE;
    get_influence_kernel<<<nblocks, NTHREAD_PER_BLOCK_INFLUENCE>>>(nspots,
                                                                   dev_endpoints,
                                                                   dev_influence,
                                                                   dev_spot_weights,
                                                                   dev_inf_volume);
    check_kernel_execution(__FILE__, __LINE__);
    retrieve_scorer<double, double4>(&(influence[0].x), dev_influence, nspots*nspots);
    retrieve_scorer<float, float>(&(inf_volume[0]), dev_inf_volume, inf_volume.size());
    retrieve_scorer<float, float>(&(spot_weights[0]), dev_spot_weights, spot_weights.size());
    // Free memory
    gpuErrchk( cudaFree(dev_endpoints) );
    gpuErrchk( cudaFree(dev_influence) );
    gpuErrchk( cudaFree(dev_inf_volume) );
    gpuErrchk( cudaFree(dev_spot_weights) );
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
