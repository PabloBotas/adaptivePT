#include "gpu_ct_to_device.cuh"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "gpu_device_globals.cuh"
#include "special_types.hpp"
#include "gpu_geometry_tools.cuh"
#include "gpu_material.cuh"
#include "gpu_errorcheck.cuh"
#include "gpu_utils.cuh"
#include "density_correction.hpp"

void gpu_ct_to_device::sendGeometries(const Volume_t& ct)
{
    gpu_ct_to_device::sendDimensions(ct);
    gpu_ct_to_device::sendDensities(ct);
    gpu_ct_to_device::sendMaterialId(ct);
}

void gpu_ct_to_device::sendDimensions(const Volume_t& ct)
//  convert external to internal geometry
{
    std::cout << "Setting CT dimensions in device..." << std::endl;
    int3   ct_n = make_int3(ct.n.x, ct.n.y, ct.n.z);
    float3 ct_d = make_float3(ct.d.x, ct.d.y, ct.d.z);

    gpuErrchk( cudaMemcpyToSymbol(ctTotalVoxN, &ct.nElements, sizeof(unsigned int), 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(ctVoxSize, &ct_d, sizeof(float3), 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(ctVox, &ct_n, sizeof(int3), 0, cudaMemcpyHostToDevice) );
}

void gpu_ct_to_device::sendDensities(const Volume_t &ct)
//  generate phantom data based on CT volume
{
    std::vector<float> densities;
    densities.reserve(ct.nElements);

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

    for(size_t i = 0; i < ct.nElements; i++)
    {
        short  val = std::max(ct.data.at(i), -1000.f);
        size_t ind = std::min(val+1000, (int)density_correction::factor.size()-1);
        densities.push_back(HU2dens(val)*density_correction::factor.at(ind));
    }

    sendVectorToTexture(ct.n.z, ct.n.y, ct.n.x, densities, dens, dens_tex);
}

void gpu_ct_to_device::sendMaterialId(const Volume_t &ct,
                                      const std::vector<int>& hu_indexes)
//  generate phantom data based on CT volume
{
    std::vector<float> materialID(ct.nElements);

    std::cout << "sendMatID: Converting HU to material ID ..." << std::endl;

    for(size_t i = 0; i < ct.nElements; i++)
    {
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
