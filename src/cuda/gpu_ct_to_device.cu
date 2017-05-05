#include "gpu_ct_to_device.cuh"

#include <iostream>
#include <string>
#include <vector>

#include "gpu_device_interaction.cuh"
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

    for(unsigned int i = 0; i < ct.n.x; i++)
    {
        for(unsigned int k=0; k < ct.n.z; k++)
        {
            for (unsigned int j=0; j < ct.n.y; j++)
            {
                int vox = getabs(i, j, k, ct.n.y, ct.n.z);
                short ctnumvox = std::max(ct.hu.at(vox), -1000.f);
                size_t ind = ctnumvox+1000;
                ind = (ind > density_correction::factor.size()-1) ? density_correction::factor.size()-1 : ind;
                densities.push_back(HU2dens(ctnumvox) * density_correction::factor.at(ind));
            }
        }
    }

    //  create a 3d array on device
    cudaExtent volumeSize = make_cudaExtent(ct.n.z, ct.n.y, ct.n.x);
    gpuErrchk( cudaMalloc3DArray(&dens, &dens_tex.channelDesc, volumeSize) );

    //  copy density data to GPU
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

//__global__ void ct_to_densities(unsigned int hu_elements, unsigned int d_elements, float* hu, float* densities, float* factor)
//{
//	const unsigned int id = blockIdx.x*blockDim.x + threadIdx.x;
//	if(id < hu_elements)
//	{
//		short hu_val = (hu[id] > -1000.f) ? hu[id] : -1000.f;
//		short ind = hu_val + 1000;
//		ind = (ind > d_elements-1) ? d_elements-1 : ind;
//		densities[id] = HU2dens(hu_val)*factor[ind];
//	}
//}


