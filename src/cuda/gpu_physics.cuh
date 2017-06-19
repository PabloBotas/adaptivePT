#ifndef __GPU_PHYSICS_CUH__
#define __GPU_PHYSICS_CUH__

__device__ void get_water_step(float& step,
                               float& step_water,
                               float& de,
                               const float max_step,
                               const float energy_in,
                               const int4& vox);

__device__ float get_residual_range (float const energy);

__device__ float get_energy_loss (float const step_water, float const energy);

__device__ float massStpRatio(const float energy, const int4& vox);

#endif