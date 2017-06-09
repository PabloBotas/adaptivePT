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
    void sendMaterialId(const Patient_Volume_t &ct);
    void sendMaterialId(const Patient_Volume_t &ct, const std::vector<int>& hu_indexes);
    void sendMassStoppingPowerRatio();
    void sendMassStoppingPowerRatio(std::vector<int>& HU_indexes);
//    void sendDensities(const Patient_Volume_t& ct,
//                            std::string densityCorrect = "/opt/utils/adaptSpotEnergies/src/phys_data/densityCorrection.dat");
}

//__global__ void ct_to_densities(unsigned int hu_elements,
//                                unsigned int d_elements,
//                                float* hu,
//                                float* densities,
//                                float* factor);

void freeCTMemory();

#endif
