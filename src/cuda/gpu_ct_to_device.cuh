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
    void sendMask(const std::vector<std::string>& mask_files,
                  const std::vector<int>& mask_importances, const CT_Dims_t& ct_dims);
    void removeMask();
    void sendMaterialId(const Volume_t &ct);
    void sendMaterialId(const Volume_t &ct, const std::vector<int>& hu_indexes);
    void sendMassStoppingPowerRatio();
    void sendMassStoppingPowerRatio(std::vector<int>& HU_indexes);
//    void sendDensities(const Volume_t& ct,
//                            std::string densityCorrect = "/opt/utils/adaptSpotEnergies/src/phys_data/densityCorrection.dat");
}

//__global__ void ct_to_densities(unsigned int hu_elements,
//                                unsigned int d_elements,
//                                float* hu,
//                                float* densities,
//                                float* factor);

void freeCTMemory();

#endif
