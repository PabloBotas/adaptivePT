#ifndef __GPU_RUN_CUH__
#define __GPU_RUN_CUH__

#include "vector2.hpp"
#include "vector4.hpp"

#include <vector>
#include <string>

void do_raytrace (const std::vector<short> spots_per_beam,
                  float4* positions_scorer,
                  float* traces_scorer = NULL);

void buffers_to_device (const std::vector<float4>& xbuffer,
                        const std::vector<float4>& vxbuffer,
                        const std::vector<short2>& ixbuffer);

void buffers_to_device (const std::vector< Vector4_t<float> >& xbuffer,
                        const std::vector< Vector4_t<float> >& vxbuffer,
                        const std::vector< Vector2_t<short> >& ixbuffer);

#endif  // GPMC_CUH
