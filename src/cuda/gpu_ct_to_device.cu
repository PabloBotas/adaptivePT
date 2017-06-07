#include "gpu_ct_to_device.cuh"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "gpu_device_globals.cuh"
#include "special_types.hpp"
#include "gpu_geometry_operations.cuh"
#include "gpu_material.cuh"
#include "gpu_errorcheck.cuh"
#include "density_correction.hpp"

void gpu_ct_to_device::sendDimensions(const Patient_Volume_t& ct)
//  convert external to internal geometry
{
    std::cout << "Setting CT dimensions in device..." << std::endl;
    int3   ct_n = make_int3(ct.n.x, ct.n.y, ct.n.z);
    float3 ct_d = make_float3(ct.d.x, ct.d.y, ct.d.z);

    gpuErrchk( cudaMemcpyToSymbol(ctTotalVoxN, &ct.nElements, sizeof(unsigned int), 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(ctVoxSize, &ct_d, sizeof(float3), 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(ctVox, &ct_n, sizeof(int3), 0, cudaMemcpyHostToDevice) );
}

void gpu_ct_to_device::sendDensities(const Patient_Volume_t &ct)
//  generate phantom data based on CT volume
{
    std::vector<float> densities;
    densities.reserve(ct.nElements);

    std::cout << "setDensities: Converting HU to densities ..." << std::endl;

//    float* gpu_densities = NULL;
//    gpuErrchk( cudaMalloc((void**) &gpu_densities, ct.nElements*sizeof(float)) );
//    float* gpu_hu = NULL;
//    gpuErrchk( cudaMalloc((void**) &gpu_hu, ct.nElements*sizeof(float)) );
//    gpuErrchk( cudaMemcpy(gpu_hu, ct.hu.data(), ct.nElements*sizeof(float), cudaMemcpyHostToDevice) );
//    float* gpu_correction = NULL;
//    gpuErrchk( cudaMalloc((void**) &gpu_correction, ct.nElements*sizeof(float)) );
//    gpuErrchk( cudaMemcpy(gpu_correction, density_correction::factor.data(), ct.nElements*sizeof(float), cudaMemcpyHostToDevice) );
//
//    int nblocks = 1 + (ct.nElements-1)/1024;
//    ct_to_densities<<<nblocks, 1024>>>(ct.nElements, density_correction::factor.size(), gpu_hu, gpu_densities, gpu_correction);
//    gpuErrchk( cudaPeekAtLastError() );
//    gpuErrchk( cudaThreadSynchronize() );
//    gpuErrchk( cudaMemcpy(&densities[0], gpu_densities, ct.nElements*sizeof(float), cudaMemcpyDeviceToHost) );

    for(size_t i = 0; i < ct.nElements; i++)
    {
        short  val = std::max(ct.hu.at(i), -1000.f);
        size_t ind = std::min(val+1000, (int)density_correction::factor.size()-1);
        densities.push_back(HU2dens(val)*density_correction::factor.at(ind));
    }

    sendVectorToTexture(ct.n.z, ct.n.y, ct.n.x, densities, dens, dens_tex);
}

void gpu_ct_to_device::sendMaterialId(const Patient_Volume_t &ct,
                                      const std::vector<int> hu_indexes)
//  generate phantom data based on CT volume
{
    std::vector<float> materialID(ct.nElements);

    std::cout << "sendMatID: Converting HU to material ID ..." << std::endl;

    for(size_t i = 0; i < ct.nElements; i++)
    {
        short val = std::max(ct.hu.at(i), -1000.f);
        materialID[i] = HU2matId(val, hu_indexes);
    }

    sendVectorToTexture(ct.n.z, ct.n.y, ct.n.x, materialID, matid, matid_tex);
}

std::vector<int> gpu_ct_to_device::sendMassStoppingPowerRatio()
{
    //  read mass stopping power ratio
    std::string file = "../src/phys_data/mass_stopping_power_ratio.dat";
    std::cout << "sendMassStoppingPowerRatio: Reading " << file << std::endl;
    std::ifstream stream(file);
    if (!stream.is_open()) {
        std::cerr << "Can't open file: " << file << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::string dummy;
    // Two dummy lines 
    std::getline(stream, line);
    std::getline(stream, line);
    // Get number of materials
    std::getline(stream, line);
    size_t const n_materials = stoi(line);

    std::vector<int> HU_starting_values(n_materials);
    std::vector<float> stp_ratios;
    size_t n_energies;
    float minimum_energy;
    float delta_energy;
    
    // Read data
    for (size_t imat = 0; imat < n_materials; imat++)
    {
        // Get number of energies per material
        std::getline(stream, line);
        std::getline(stream, line);
        std::istringstream ss(line);
        ss >> dummy >> HU_starting_values.at(imat) >>
              dummy >> n_energies >> minimum_energy >> delta_energy;
        std::getline(stream, line);

        if (imat == 0)
        {
            stp_ratios.resize(n_energies*n_materials);
        }

        for (size_t i = 0; i < n_energies; i++)
        {
            std::getline(stream, line);
            std::istringstream ss(line);
            ss >> dummy >> stp_ratios.at(i + imat*n_energies);
        }
    }

    //  transfer to GPU
    minimum_energy *= MeV2eV;
    delta_energy *= MeV2eV;
    cudaMemcpyToSymbol(stp_ratio_min_e, &minimum_energy, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(stp_ratio_delta_e, &delta_energy, sizeof(float), 0, cudaMemcpyHostToDevice);

    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaMallocArray(&stp_ratio_array, &channelDesc, n_energies, n_materials);
    cudaMemcpyToArray(stp_ratio_array, 0,0, &stp_ratios[0], sizeof(float)*stp_ratios.size(), cudaMemcpyHostToDevice);
    stp_ratio_tex.filterMode = cudaFilterModeLinear;
    cudaBindTextureToArray(stp_ratio_tex, stp_ratio_array, channelDesc);

    return HU_starting_values;
}

void sendVectorToTexture (size_t w, size_t h, size_t d,
                          std::vector<float> host_vec,
                          cudaArray* array,
                          texture<float, 3, cudaReadModeElementType>& tex)
{
    //  create a 3d array on device
    cudaExtent extent = make_cudaExtent(w, h, d);
    gpuErrchk( cudaMalloc3DArray(&array, &tex.channelDesc, extent) );

    // copy data to GPU
    cudaMemcpy3DParms pars = {0};
    pars.srcPtr   = make_cudaPitchedPtr((void *)host_vec.data(),
                                        extent.width*sizeof(float),
                                        extent.width, extent.height);
    pars.dstArray = array;
    pars.extent   = extent;
    pars.kind     = cudaMemcpyHostToDevice;
    gpuErrchk( cudaMemcpy3D(&pars) );
    // Bind device array to texture
    tex.normalized = false;
    tex.filterMode = cudaFilterModePoint;
    gpuErrchk( cudaBindTextureToArray(tex, array, tex.channelDesc) );
}


//__global__ void ct_to_densities(unsigned int hu_elements, unsigned int d_elements, float* hu, float* densities, float* factor)
//{
//    const unsigned int id = blockIdx.x*blockDim.x + threadIdx.x;
//    if(id < hu_elements)
//    {
//        short hu_val = (hu[id] > -1000.f) ? hu[id] : -1000.f;
//        short ind = hu_val + 1000;
//        ind = (ind > d_elements-1) ? d_elements-1 : ind;
//        densities[id] = HU2dens(hu_val)*factor[ind];
//    }
//}
