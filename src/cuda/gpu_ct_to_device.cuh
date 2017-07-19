#ifndef __GPU_CT_TO_DEVICE_CUH__
#define __GPU_CT_TO_DEVICE_CUH__

#include "volume.hpp"
#include "special_types.hpp"

#include <string>
#include <vector>

namespace gpu_ct_to_device
{
    void sendGeometries(const Volume_t& ct);
    void sendDimensions(const Volume_t& ct);
    void sendDensities(const Volume_t& ct);
    void sendMaterialId(const Volume_t &ct);
    void sendMaterialId(const Volume_t &ct, const std::vector<int>& hu_indexes);
    void sendMassStoppingPowerRatio();
    void sendMassStoppingPowerRatio(std::vector<int>& HU_indexes);
//    void sendDensities(const Volume_t& ct,
//                            std::string densityCorrect = "/opt/utils/adaptSpotEnergies/src/phys_data/densityCorrection.dat");
}

//__global__ void ct_to_densities(unsigned int hu_elements,
//                                unsigned int d_elements,
//                                double* hu,
//                                double* densities,
//                                double* factor);

void freeCTMemory();

#endif
