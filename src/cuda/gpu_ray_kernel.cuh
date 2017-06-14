#ifndef __GPU_RAY_KERNEL_CUH__
#define __GPU_RAY_KERNEL_CUH__

__global__ void raytrace_plan_kernel(const int num,
                                     const short *spots_per_beam,
                                     float4 *pos_scorer,
                                     float* traces = NULL);

__device__ unsigned int get_endpoints_index(const short beam_id,
                                            const short spot_id,
                                            const short* spots_per_beam);

__device__ void getWaterStep(float& step,
                             float& step_water,
                             const float max_step,
                             const float energy,
                             const float avail_wepl,
                             const int4& vox);

__device__ float massStpRatio(const float energy, const int4& vox);

#endif
