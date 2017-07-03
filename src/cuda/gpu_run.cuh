#ifndef __GPU_RUN_CUH__
#define __GPU_RUN_CUH__

#include "vector2.hpp"
#include "vector4.hpp"

#include <vector>
#include <string>

template<class T>
void do_raytrace (const std::vector<short>& spots_per_beam,
                  float4* positions_scorer,
                  T* traces_scorer,
                  const Array4<float>& orig_endpoints);

void buffers_to_device (const std::vector<float4>& xbuffer,
                        const std::vector<float4>& vxbuffer,
                        const std::vector<short2>& ixbuffer,
                        const bool alloc = true);

#endif  // GPMC_CUH
