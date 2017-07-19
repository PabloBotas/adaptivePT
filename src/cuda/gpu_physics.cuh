#ifndef __GPU_PHYSICS_CUH__
#define __GPU_PHYSICS_CUH__

__device__ void get_water_step(double& step,
                               double& step_water,
                               double& de,
                               const double max_step,
                               const double energy_in,
                               const int4& vox);

__device__ double get_residual_range (double const energy);

__device__ double get_energy_loss (double const step_water, double const energy);

__device__ double massStpRatio(const double energy, const int4& vox);

#endif