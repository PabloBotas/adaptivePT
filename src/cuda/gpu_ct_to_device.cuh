#ifndef __GPU_CT_TO_DEVICE_CUH__
#define __GPU_CT_TO_DEVICE_CUH__

#include "volume.hpp"
#include "special_types.hpp"

#include <string>
#include <vector>

namespace gpu_ct_to_device
{
    void sendGeometries(const Patient_Volume_t& ct);
    void sendDimensions(const Patient_Volume_t& ct);
    void sendDensities(const Patient_Volume_t& ct);
    void sendMaterialId(const Patient_Volume_t &ct, const std::vector<int> hu_indexes);
    std::vector<int> sendMassStoppingPowerRatio();
//    void sendDensities(const Patient_Volume_t& ct,
//                            std::string densityCorrect = "/opt/utils/adaptSpotEnergies/src/phys_data/densityCorrection.dat");
}

void sendVectorToTexture(size_t w, size_t h, size_t d,
                         std::vector<float> host_vec,
                         cudaArray* array,
                         texture<float, 3, cudaReadModeElementType>& tex);

//__global__ void ct_to_densities(unsigned int hu_elements,
//                                unsigned int d_elements,
//                                float* hu,
//                                float* densities,
//                                float* factor);

void freeCTMemory();

#endif
