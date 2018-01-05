#include "gpu_influence_launchers.cuh"

#include "gpu_device_globals.cuh"
#include "gpu_influence_kernel.cuh"
#include "gpu_utils.cuh"
#include "patient_parameters.hpp"
#include "special_types.hpp"
#include "tramp.hpp"
#include "vector4.hpp"
#include "volume.hpp"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <string>
#include <vector>

void influence_from_beam_model_launcher(std::string outfile,
                                        Array4<float>& influence,
                                        const std::vector<float>& new_energies,
                                        const Volume_metadata_t& ct_metadata,
                                        const Patient_Parameters_t& patient_parameters,
                                        const uint nspots, const uint n_probing_positions)
{
    std::vector<float> inf_cube(ct_metadata.nElements);
    std::vector<float> spot_weights;
    for (uint i = 0; i < patient_parameters.tramp_files.size(); ++i) {
        Tramp_t tramp(patient_parameters.tramp_files.at(i), patient_parameters.machine);
        std::vector<float> w = tramp.get_weights();
        spot_weights.reserve(spot_weights.size() + tramp.nspots);
        spot_weights.insert(spot_weights.end(), w.begin(), w.end());
    }

    // Create scorer array
    float4 *dev_influence = NULL;
    array_to_device<float4, Vector4_t<float> >(dev_influence, influence.data(), influence.size());
    // Create scorer array
    float *dev_inf_volume = NULL;
    allocate_scorer<float>(dev_inf_volume, inf_cube.size());
    // Create weights array
    float *dev_spot_weights = NULL;
    array_to_device<float>(dev_spot_weights, &spot_weights[0], spot_weights.size());
    // Create new energies array
    float *dev_new_energies = NULL;
    if (new_energies.size() > 0)
        array_to_device<float>(dev_new_energies, &new_energies[0], new_energies.size());
    
    // Launch influence kernel
    uint matrix_elements = nspots*n_probing_positions;
    std::cout << "Calculating influence matrix with " << matrix_elements << " elements (";
    std::cout << (float)matrix_elements*sizeof(float)/1E9 << " Gb) ..." << std::endl;

    uint nblocks = 1 + (matrix_elements-1)/NTHREAD_PER_BLOCK_INFLUENCE;
    get_influence_kernel<<<nblocks, NTHREAD_PER_BLOCK_INFLUENCE>>>(nspots, n_probing_positions,
                                                                   dev_influence,
                                                                   dev_spot_weights,
                                                                   dev_inf_volume,
                                                                   dev_new_energies);
    check_kernel_execution(__FILE__, __LINE__);
    retrieve_scorer<float, float4>(&(influence[0].x), dev_influence, influence.size());
    retrieve_scorer<float, float>(&(inf_cube[0]), dev_inf_volume, inf_cube.size());
    // Free memory
    gpuErrchk( cudaFree(dev_influence) );
    gpuErrchk( cudaFree(dev_inf_volume) );
    gpuErrchk( cudaFree(dev_spot_weights) );
    gpuErrchk( cudaFree(dev_new_energies) );

    std::cout << "Writting Dij from beam model to: ";
    std::cout << outfile << std::endl;
    std::ofstream fout(outfile, std::ios::out | std::ios::binary);
    for (uint i = 0; i < influence.size(); ++i)
        fout.write((char*)&influence[i].w, sizeof(float));

#ifdef __DEBUG_INFLUENCE_MATRICES__
    bool ct_case = new_energies.size() == 0;
    std::string out_case = ct_case ? "CT" : "CBCT";
    std::cout << "Writting influence_cube_"+out_case+".dat ..." << std::endl;
    std::outdir = utils::get_full_parent_path(outfile);
    std::string file = outdir+"/influence_cube_"+out_case+".dat";
    std::ofstream fout2(file, std::ios::out | std::ios::binary);
    fout2.write((char*)inf_cube.data(), inf_cube.size()*sizeof(float));

    std::cout << "Writting vox_endpoints_"+out_case+".dat ..." << std::endl;
    file = outdir+"/vox_endpoints_"+out_case+".dat";
    std::ofstream fout3(file, std::ios::out | std::ios::binary);
    for (uint i = 0; i < nspots; ++i) {
        unsigned int vox_x = floor(influence.at(i).x/ct_metadata.d.x);
        unsigned int vox_y = floor(influence.at(i).y/ct_metadata.d.y);
        unsigned int vox_z = floor(influence.at(i).z/ct_metadata.d.z);
        int vox_w;
        // Check if in CT grid
        if (vox_x >= ct_metadata.n.x ||
            vox_y >= ct_metadata.n.y ||
            vox_z >= ct_metadata.n.z)
            vox_w = -1;
        else
            vox_w = vox_z + vox_y*ct_metadata.n.z + vox_x*ct_metadata.n.z*ct_metadata.n.y;

        fout3.write((char*)&vox_w, sizeof(int));
    }
    // for (size_t i = 0; i < 10; i++) {
    //     std::cerr << influence[i].x << " " << influence[i].y << " " << influence[i].z << " " << influence[i].w << "\n";
    // }
#endif
}