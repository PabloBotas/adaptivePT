#ifndef __PHYSICS_TO_DEVICE_CUH__
#define __PHYSICS_TO_DEVICE_CUH__

#include <vector>

namespace gpu_physics_to_device
{
    void sendMassStoppingPowerRatio();
    void sendMassStoppingPowerRatio(std::vector<int>& HU_starting_values);
    void sendWaterRestrictedSPower();
}
void freePhysicsMemory();

#endif