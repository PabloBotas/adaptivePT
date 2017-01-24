#ifndef __GPU_RUN_CUH__
#define __GPU_RUN_CUH__

#include "patient_parameters.hpp"
#include "special_types.hpp"

#include <vector>
#include <string>

void init_rays(const Patient_Parameters_t& pat,
               const unsigned int ibeam,
               std::vector<float4>& xbuffer,
               std::vector<float4>& vxbuffer);

void calculateRays(std::vector<float4>& xbuffer,
                   std::vector<float4>& vxbuffer,
                   const BeamAngles_t& ang,
                   const float3& ct_offsets);

unsigned int rays_to_device(std::vector<float4>& xbuffer,
                            std::vector<float4>& vxbuffer,
                            const float2& angles,
                            const float3& ct_offsets);

__global__ void rays_to_device_kernel(int num, float2 angles, float3 ct_offsets);

void clearScorer(void *s, size_t sz);

// void outputData(void *src, const size_t size, string outputfilename, const char *mode);
// void outputScorerResults(size_t beamNumber);
void freeMemory();

#endif  // GPMC_CUH
