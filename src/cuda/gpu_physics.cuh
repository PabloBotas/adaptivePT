#ifndef __GPU_PHYSICS_CUH__
#define __GPU_PHYSICS_CUH__

#include "gpu_ray_class.cuh"

__device__ void get_step(double& step,
                         double& step_water,
                         double& de,
                         const double max_step,
                         const double energy_in,
                         const int4& vox);

// Q50 blur
__device__ void get_q50_blur_step(double& step,
                                  double& step_water,
                                  double& de,
                                  const double max_step,
                                  const Ray& ray,
                                  const int4& vox);
__device__ void q50_density_blur(double& density, int4& q50_vox,
                                 const Ray& ray, const int4& vox);
// Average blur
__device__ void get_average_blur_step(double& step,
                                      double& step_water,
                                      double& de,
                                      const double max_step,
                                      const Ray& ray,
                                      const int4& vox);
__device__ void average_blur (double& density, double& mass_stp_ratio,
                              const Ray& ray, const int4& vox);

__device__ void step_blur_get_voxels(int* voxels, const Ray& ray, const int4& vox,
                                     const int n, const double scale);
__device__ double get_residual_range (double const energy);

__device__ double get_energy_loss (double const step_water, double const energy);

__device__ double massStpRatio(const double energy, const int4& vox);

__device__ double get_dose_at (double const energy, double const wepl_r, double const wepl_d);

#endif