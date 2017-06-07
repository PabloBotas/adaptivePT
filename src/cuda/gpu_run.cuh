#ifndef __GPU_RUN_CUH__
#define __GPU_RUN_CUH__

#include "special_types.hpp"

#include <vector>
#include <string>

void calculateRays(const std::vector<short> spots_per_beam,
                   float4* positions_scorer,
                   float4* directions_scorer,
                   short2* metadata_scorer,
                   float* traces_scorer = NULL);

void buffers_to_device(const std::vector<float4>& xbuffer,
                       const std::vector<float4>& vxbuffer,
                       const std::vector<short2>& ixbuffer);

#endif  // GPMC_CUH
