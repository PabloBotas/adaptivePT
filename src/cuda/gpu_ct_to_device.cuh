#ifndef __GPU_CT_TO_DEVICE_CUH__
#define __GPU_CT_TO_DEVICE_CUH__

#include "volume.hpp"
#include "special_types.hpp"

#include <string>

namespace gpu_ct_to_device
{
    void setDimensions(const Patient_Volume_t& ct);
    void setDensities(const Patient_Volume_t& ct,
                      std::string densityCorrect = "src/phys_data/densityCorrection.dat");
}

#endif