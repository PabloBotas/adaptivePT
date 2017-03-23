#include "gpu_ct_to_device.cuh"

#include <iostream>
#include <string>
#include <vector>

#include "gpu_device_interaction.cuh"
#include "special_types.hpp"
#include "gpu_geometry_operations.cuh"
#include "gpu_material.cuh"
#include "gpu_errorcheck.cuh"

void gpu_ct_to_device::setDimensions(const Patient_Volume_t& ct)
//  convert external to internal geometry
{
    std::cout << "Setting CT dimensions in device..." << std::endl;
    int3   ct_n = make_int3(ct.n.x, ct.n.y, ct.n.z);
    float3 ct_d = make_float3(ct.d.x, ct.d.y, ct.d.z);

    // Swap CT grid X, Z axis
    std::swap(ct_d.x, ct_d.z);
    std::swap(ct_n.x, ct_n.z);

    gpuErrchk( cudaMemcpyToSymbol(ctTotalVoxN, &ct.nElements, sizeof(unsigned int), 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(ctVoxSize, &ct_d, sizeof(float3), 0, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpyToSymbol(ctVox, &ct_n, sizeof(int3), 0, cudaMemcpyHostToDevice) );
}


void gpu_ct_to_device::setDensities(const Patient_Volume_t &ct,
                                    std::string densityCorrect)
//  generate phantom data based on CT volume
{
    std::vector<float> densities;
    densities.reserve(ct.nElements);
    std::vector<float> dcfactor = readDensityCorrect(densityCorrect);

    std::cout << "densities_to_device: Converting CT to material phantom..." << std::endl;
    for(unsigned int i = 0; i < ct.n.x; i++)
    {
        for(unsigned int k=0; k < ct.n.z; k++)
        {
            for (unsigned int j=0; j < ct.n.y; j++)
            {
                int vox = getabs(i, j, k, ct.n.y, ct.n.z);
                short ctnumvox = std::max(ct.hu.at(vox), -1000.f);
                size_t ind = ctnumvox+1000;
                ind = (ind > dcfactor.size()-1) ? dcfactor.size()-1 : ind;
                densities.push_back(HU2dens(ctnumvox) * dcfactor.at(ind));              
            }
        }
    }

    //  create a 3d array on device
    cudaExtent volumeSize = make_cudaExtent(ct.n.z, ct.n.y, ct.n.x);
    gpuErrchk( cudaMalloc3DArray(&dens, &dens_tex.channelDesc, volumeSize) );

    //  copy denstiy data to GPU
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr = make_cudaPitchedPtr((void*)&densities[0], volumeSize.width*sizeof(float),
                                            volumeSize.width, volumeSize.height);
    copyParams.dstArray = dens;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;
    gpuErrchk( cudaMemcpy3D(&copyParams) );
    //  copy data from host to device
    dens_tex.normalized = false;
    dens_tex.filterMode = cudaFilterModePoint;
    gpuErrchk( cudaBindTextureToArray(dens_tex, dens, dens_tex.channelDesc) );
    //  bind to texture memory
}

