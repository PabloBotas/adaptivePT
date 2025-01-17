#include "gpu_ct_to_device.cuh"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "gpu_device_globals.cuh"
#include "special_types.hpp"
#include "gpu_geometry_tools.cuh"
#include "gpu_material.cuh"
#include "gpu_errorcheck.cuh"
#include "gpu_utils.cuh"
#include "density_correction.hpp"
#include "utils.hpp"

void gpu_ct_to_device::sendGeometries(const Volume_t& ct)
{
    gpu_ct_to_device::sendDimensions(ct);
    gpu_ct_to_device::sendDensities(ct);
    gpu_ct_to_device::sendMaterialId(ct);
}

void gpu_ct_to_device::sendDimensions(const Volume_t& ct)
//  convert external to internal geometry
{
    std::cout << "sendDimensions: Setting volume dimensions in device..." << std::endl;
    int3   ct_n = make_int3(ct.n.x, ct.n.y, ct.n.z);
    float3 ct_d = make_float3(ct.d.x, ct.d.y, ct.d.z);

    gpuErrchk( cudaMemcpyToSymbol(ctTotalVoxN, &ct.nElements, sizeof(unsigned int), 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(ctVoxSize, &ct_d, sizeof(float3), 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(ctVox, &ct_n, sizeof(int3), 0, cudaMemcpyHostToDevice) );
}

void gpu_ct_to_device::sendDensities(const Volume_t &ct)
//  generate phantom data based on CT volume
{
    std::vector<float> densities(ct.nElements);

    std::cout << "setDensities: Converting HU to densities ..." << std::endl;

//    float* gpu_densities = NULL;
//    gpuErrchk( cudaMalloc((void**) &gpu_densities, ct.nElements*sizeof(float)) );
//    float* gpu_hu = NULL;
//    gpuErrchk( cudaMalloc((void**) &gpu_hu, ct.nElements*sizeof(float)) );
//    gpuErrchk( cudaMemcpy(gpu_hu, ct.data.data(), ct.nElements*sizeof(float), cudaMemcpyHostToDevice) );
//    float* gpu_correction = NULL;
//    gpuErrchk( cudaMalloc((void**) &gpu_correction, ct.nElements*sizeof(float)) );
//    gpuErrchk( cudaMemcpy(gpu_correction, density_correction::factor.data(), ct.nElements*sizeof(float), cudaMemcpyHostToDevice) );
//
//    int nblocks = 1 + (ct.nElements-1)/1024;
//    ct_to_densities<<<nblocks, 1024>>>(ct.nElements, density_correction::factor.size(), gpu_hu, gpu_densities, gpu_correction);
//    gpuErrchk( cudaPeekAtLastError() );
//    gpuErrchk( cudaThreadSynchronize() );
//    gpuErrchk( cudaMemcpy(&densities[0], gpu_densities, ct.nElements*sizeof(float), cudaMemcpyDeviceToHost) );

    #pragma omp parallel for
    for (size_t i = 0; i < ct.nElements; i++) {
        short  val = std::max(ct.data.at(i), -1000.f);
        size_t ind = std::min(val+1000, (int)density_correction::factor.size()-1);
        densities.at(i) = HU2dens(val)*density_correction::factor.at(ind);
    }

    sendVectorToTexture(ct.n.z, ct.n.y, ct.n.x, densities, dens, dens_tex);
}

void gpu_ct_to_device::sendMaterialId(const Volume_t &ct,
                                      const std::vector<int>& hu_indexes)
//  generate phantom data based on CT volume
{
    std::vector<float> materialID(ct.nElements);

    std::cout << "sendMatID: Converting HU to material ID ..." << std::endl;

    #pragma omp parallel for
    for (size_t i = 0; i < ct.nElements; i++) {
        short val = std::max(ct.data.at(i), -1000.f);
        materialID[i] = HU2matId(val, hu_indexes);
    }

    sendVectorToTexture(ct.n.z, ct.n.y, ct.n.x, materialID, matid, matid_tex);
}

void gpu_ct_to_device::sendMaterialId(const Volume_t &ct)
{
    std::vector<int> hu_indexes = {
        -1000, -950, -120, -83, -53, -23, 7,
        18, 80, 120, 200, 300, 400, 500, 600,
        700, 800, 900, 1000, 1100, 1200, 1300,
        1400, 1500, 2995, 0};

    gpu_ct_to_device::sendMaterialId(ct, hu_indexes);
}


void freeCTMemory()
{
    cudaFreeArray(dens);
    cudaUnbindTexture(dens_tex);
    cudaFreeArray(matid);
    cudaUnbindTexture(matid_tex);
}


void gpu_ct_to_device::sendMask(const std::vector<std::string>& mask_files,
                                const std::vector<int>& mask_importances,
                                const CT_Dims_t& ct_dims)
//  generate phantom data based on CT volume
{
    bool ifmasks = mask_files.size() ? true : false;
    gpuErrchk( cudaMemcpyToSymbol(masking_vf, &ifmasks, sizeof(bool), 0, cudaMemcpyHostToDevice) );

    if (ifmasks) {
        const float mask_threshold = 0.5;
        Volume_t vol = utils::read_masks (mask_files, mask_threshold, mask_importances);
        if (vol.nElements != ct_dims.total) {
            std::cerr << "ERROR! The number of elements found in the VF mask differs from ";
            std::cerr << "the number of elements in the CT!!" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::vector<int> vec(vol.nElements);
        std::transform (vol.data.begin(), vol.data.end(), vec.begin(),
            [mask_threshold](const float& d) {return d > mask_threshold ? int(d+0.5) : 0;});

        std::ofstream ofs("vf_mask_map.dat", std::ios::binary);
        ofs.write((char*)vec.data(), vec.size()*sizeof(int));
        ofs.close();

        std::cout << "sendMask: Sending VF mask to device ..." << std::endl;
        sendVectorToTexture(ct_dims.n.z, ct_dims.n.y, ct_dims.n.x, vec, vf_mask, vf_mask_tex);
    }
}


void gpu_ct_to_device::removeMask()
//  generate phantom data based on CT volume
{
    bool ifmasks = false;
    gpuErrchk( cudaMemcpyToSymbol(masking_vf, &ifmasks, sizeof(bool), 0, cudaMemcpyHostToDevice) );
    cudaFreeArray(vf_mask);
    cudaUnbindTexture(vf_mask_tex);
}


//__global__ void ct_to_densities(unsigned int hu_elements, unsigned int d_elements, float* data, float* densities, float* factor)
//{
//    const unsigned int id = blockIdx.x*blockDim.x + threadIdx.x;
//    if(id < hu_elements)
//    {
//        short hu_val = (data[id] > -1000.f) ? data[id] : -1000.f;
//        short ind = hu_val + 1000;
//        ind = (ind > d_elements-1) ? d_elements-1 : ind;
//        densities[id] = HU2dens(hu_val)*factor[ind];
//    }
//}
