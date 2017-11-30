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
                            const Volume_t& ct,
                            Array4<float>& endpoints,
                            Array4<float>& initpos_xbuffer_dbg,
                            Array4<float>& initpos,
                            const Parser& parser,
                            Array4<float>& influence)
{
    // Set geometry in GPU
    gpu_ct_to_device::sendGeometries(ct);

    // Create host buffers and initialize rays
    std::vector<float4> xbuffer;
    std::vector<float4> vxbuffer;
    std::vector<short2> ixbuffer;
    create_virtual_source_buffers (pat, xbuffer, vxbuffer, ixbuffer);
    buffers_to_device (xbuffer, vxbuffer, ixbuffer, true);
    virtual_src_to_treatment_plane (xbuffer.size(), pat.angles,
                                    make_float3(pat.ct.offset.x, pat.ct.offset.y, pat.ct.offset.z));

    gpuErrchk( cudaMemcpyFromSymbol(&(initpos[0].x), xdata, sizeof(float4)*xbuffer.size(), 0, cudaMemcpyDeviceToHost) );
    // Copy buffer with initial positions and wepl
    for (size_t i = 0; i < xbuffer.size(); i++) {
        initpos_xbuffer_dbg.at(i).x = xbuffer[i].x;
        initpos_xbuffer_dbg.at(i).y = xbuffer[i].y;
        initpos_xbuffer_dbg.at(i).z = xbuffer[i].z;
        initpos_xbuffer_dbg.at(i).w = xbuffer[i].w; // wepl
    }

    gpu_raytrace (pat, endpoints, parser.output_ct_traces);

    // INFLUENCE
    std::vector<float> inf_cube(ct.nElements);
    std::vector<float> spot_weights;
    for (size_t i = 0; i < pat.tramp_files.size(); ++i) {
        Tramp_t tramp(pat.tramp_files.at(i), pat.machine);
        std::vector<float> w = tramp.get_weights();
        spot_weights.reserve(spot_weights.size() + tramp.nspots);
        spot_weights.insert(spot_weights.end(), w.begin(), w.end());
    }
    uint nprobes = uint(influence.size()/pat.total_spots + 0.5);
    gpu_calculate_influence (pat.total_spots, nprobes, influence, spot_weights, inf_cube);

#ifdef __INFLUENCE_MATRICES__
    std::cout << "Writting influence_matrix_CT.dat ..." << std::endl;
    std::string outputdir = parser.output_opt4D_files.empty() ? parser.out_dir :
                                                                parser.output_opt4D_files;
    std::ofstream fout(outputdir+"/influence_matrix_CT.dat", std::ios::out | std::ios::binary);
    for (size_t i = 0; i < influence.size(); ++i)
        fout.write((char*)&influence[i].w, sizeof(float));

    std::cout << "Writting influence_cube_CT.dat ..." << std::endl;
    std::ofstream fout2(outputdir+"/influence_cube_CT.dat", std::ios::out | std::ios::binary);
    fout2.write((char*)inf_cube.data(), inf_cube.size()*sizeof(float));

    std::cout << "Writting vox_endpoints_CT.dat ..." << std::endl;
    std::ofstream fout3(outputdir+"/vox_endpoints_CT.dat", std::ios::out | std::ios::binary);
    for (size_t i = 0; i < pat.total_spots; ++i) {
        unsigned int vox_x = floor(influence.at(i).x/ct.d.x);
        unsigned int vox_y = floor(influence.at(i).y/ct.d.y);
        unsigned int vox_z = floor(influence.at(i).z/ct.d.z);
        int vox_w;
        // Check if in CT grid
        if (vox_x >= ct.n.x ||
            vox_y >= ct.n.y ||
            vox_z >= ct.n.z)
            vox_w = -1;
        else
            vox_w = vox_z + vox_y*ct.n.z + vox_x*ct.n.z*ct.n.y;

        fout3.write((char*)&vox_w, sizeof(int));
    }
#endif
    // for (size_t i = 0; i < 10; i++)
    // {
    //     std::cerr << influence[i].x << " " << influence[i].y << " " << influence[i].z << " " << influence[i].w << "\n";
    // }
    freeCTMemory();
}

void gpu_raytrace_warped (const Patient_Parameters_t &pat,
                          const Volume_t &ct,
                          const Array4<float>& orig_endpoints,
                          const Array4<float>& init_pos,
                          Array4<float>& endpoints,
                          const Parser& parser,
                          Array4<float>& influence)
{
    // Set geometry in GPU
    gpu_ct_to_device::sendGeometries(ct);

    // Create host buffers and initialize rays
    std::vector<float4> xbuffer;
    std::vector<float4> vxbuffer;
    std::vector<short2> ixbuffer;
    create_treatment_plane_buffers (pat, orig_endpoints, init_pos,
                                    xbuffer, vxbuffer, ixbuffer);
    buffers_to_device (xbuffer, vxbuffer, ixbuffer, false);
    correct_offsets (xbuffer.size(), 
                     make_float3(pat.ct.offset.x, pat.ct.offset.y, pat.ct.offset.z),
                     make_float3(pat.original_ct.offset.x, pat.original_ct.offset.y, pat.original_ct.offset.z));
    Array4<float> off_endpoints = offset_endpoints (orig_endpoints, 
                                                     make_float3(pat.ct.offset.x, pat.ct.offset.y, pat.ct.offset.z),
                                                     make_float3(pat.original_ct.offset.x, pat.original_ct.offset.y, pat.original_ct.offset.z));

    gpu_raytrace (pat, endpoints, parser.output_cbct_traces, off_endpoints);

    std::vector<float> inf_cube(ct.nElements);
    std::vector<float> spot_weights;
    for (size_t i = 0; i < pat.tramp_files.size(); ++i) {
        Tramp_t tramp(pat.tramp_files.at(i), pat.machine);
        std::vector<float> w = tramp.get_weights();
        spot_weights.reserve(spot_weights.size() + tramp.nspots);
        spot_weights.insert(spot_weights.end(), w.begin(), w.end());
    }
    std::vector<float> new_energies(endpoints.size());
    for (size_t i = 0; i < new_energies.size(); ++i) {
        new_energies.at(i) = vxbuffer.at(i).w + endpoints.at(i).w;
    }

    uint nprobes = uint(influence.size()/pat.total_spots + 0.5);
    gpu_calculate_influence (pat.total_spots, nprobes, influence,
                             spot_weights, inf_cube, new_energies);
#ifdef __INFLUENCE_MATRICES__
    std::cout << "Writting influence_CBCT.dat ..." << std::endl;
    std::string outputdir = parser.output_opt4D_files.empty() ? parser.out_dir :
                                                            parser.output_opt4D_files;
    std::ofstream fout(outputdir+"/influence_matrix_CBCT.dat", std::ios::out | std::ios::binary);
    for (size_t i = 0; i < influence.size(); ++i)
        fout.write((char*)&influence[i].w, sizeof(float));
    
    std::cout << "Writting influence_cube_CBCT.dat ..." << std::endl;
    std::ofstream fout2(outputdir+"/influence_cube_CBCT.dat", std::ios::out | std::ios::binary);
    fout2.write((char*)inf_cube.data(), inf_cube.size()*sizeof(float));

    std::cout << "Writting vox_endpoints_CBCT.dat ..." << std::endl;
    std::ofstream fout3(outputdir+"/vox_endpoints_CBCT.dat", std::ios::out | std::ios::binary);
    for (size_t i = 0; i < pat.total_spots; ++i) {
        unsigned int vox_x = floor(influence.at(i).x/ct.d.x);
        unsigned int vox_y = floor(influence.at(i).y/ct.d.y);
        unsigned int vox_z = floor(influence.at(i).z/ct.d.z);
        int vox_w;
        // Check if in CT grid
        if (vox_x >= ct.n.x ||
            vox_y >= ct.n.y ||
            vox_z >= ct.n.z)
            vox_w = -1;
        else
            vox_w = vox_z + vox_y*ct.n.z + vox_x*ct.n.z*ct.n.y;

        fout3.write((char*)&vox_w, sizeof(int));
    }
    // for (size_t i = 0; i < influence.size(); i+=6421)
    // {
    //     std::cout << influence[i].x << " " << influence[i].y << " " << influence[i].z << " " << influence[i].w << "\n";
    // }
#endif
    freeCTMemory();
}

void gpu_raytrace (const Patient_Parameters_t& pat,
                   Array4<float>& endpoints,
                   std::string traces_file,
                   const Array4<float>& orig_endpoints)
{
    // Create scorer array
    float4* pos_scorer = NULL;
    allocate_scorer<float4>(pos_scorer, pat.total_spots);

    // Calculate rays
    if (traces_file.empty()) {
        do_raytrace(pat.spots_per_field, pos_scorer, NULL, orig_endpoints);
    } else {
        float* traces_scorer = NULL;
        allocate_scorer<float>(traces_scorer, pat.ct.total);
        do_raytrace(pat.spots_per_field, pos_scorer, traces_scorer, orig_endpoints);
        bool long_data = false;
        Volume_t traces(pat.ct, long_data);
        retrieve_scorer<float, float>(&traces.data[0], traces_scorer, traces.nElements);
        // traces.output("output_volume.mha");
        traces.output(traces_file);
        gpuErrchk( cudaFree(traces_scorer) );
    }

    retrieve_scorer<float, float4>(&(endpoints[0].x), pos_scorer, pat.total_spots);
    // Free memory
    gpuErrchk( cudaFree(pos_scorer) );
}

void gpu_calculate_influence (const uint& nspots,
                              const uint& nprobes,
                              Array4<float>& influence,
                              std::vector<float>& spot_weights,
                              std::vector<float>& inf_volume,
                              const std::vector<float>& new_energies)
{
    // Create scorer array
    float4 *dev_influence = NULL;
    array_to_device<float4, Vector4_t<float> >(dev_influence, influence.data(), influence.size());
    // Create scorer array
    float *dev_inf_volume = NULL;
    allocate_scorer<float>(dev_inf_volume, inf_volume.size());
    // Create weights array
    float *dev_spot_weights = NULL;
    array_to_device<float>(dev_spot_weights, &spot_weights[0], spot_weights.size());
    // Create new energies array
    float *dev_new_energies = NULL;
    array_to_device<float>(dev_new_energies, &new_energies[0], new_energies.size());
    
    // Launch influence kernel
    std::cout << "Calculating influence matrix with " << nspots*nprobes << " elements ..." << std::endl;

    uint nblocks = 1 + (nspots*nprobes-1)/NTHREAD_PER_BLOCK_INFLUENCE;
    get_influence_kernel<<<nblocks, NTHREAD_PER_BLOCK_INFLUENCE>>>(nspots, nprobes,
                                                                   dev_influence,
                                                                   dev_spot_weights,
                                                                   dev_inf_volume,
                                                                   dev_new_energies);
    check_kernel_execution(__FILE__, __LINE__);
    retrieve_scorer<float, float4>(&(influence[0].x), dev_influence, influence.size());
    retrieve_scorer<float, float>(&(inf_volume[0]), dev_inf_volume, inf_volume.size());
    retrieve_scorer<float, float>(&(spot_weights[0]), dev_spot_weights, spot_weights.size());
    // Free memory
    gpuErrchk( cudaFree(dev_influence) );
    gpuErrchk( cudaFree(dev_inf_volume) );
    gpuErrchk( cudaFree(dev_spot_weights) );
    gpuErrchk( cudaFree(dev_new_energies) );
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
