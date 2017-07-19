#ifndef __GPU_RUN_CUH__
#define __GPU_RUN_CUH__

#include "vector2.hpp"
#include "vector4.hpp"

#include <vector>
#include <string>

void do_raytrace (const std::vector<short>& spots_per_beam,
                  double4* positions_scorer,
                  float* traces_scorer,
                  const Array4<double>& orig_endpoints);

void buffers_to_device (const std::vector<double4>& xbuffer,
                        const std::vector<double4>& vxbuffer,
                        const std::vector<short2>& ixbuffer,
                        const bool alloc = true);

#endif  // GPMC_CUH
