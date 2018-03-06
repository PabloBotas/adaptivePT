#ifndef __GPU_PHYSICS_CUH__
#define __GPU_PHYSICS_CUH__

#include "gpu_ray_class.cuh"

__device__ void get_step(float& step,
                         float& step_water,
                         float& de,
                         const float max_step,
                         const Ray& ray,
                         const int4& vox);

__device__ void get_step(float& step,
                         float& step_water,
                         float& de,
                         const float max_step,
                         const float energy_in,
                         const int4& vox);

// Q50 blur
__device__ void get_q50_blur_step(float& step,
                                  float& step_water,
                                  float& de,
                                  const float max_step,
                                  const Ray& ray,
                                  const int4& vox);
__device__ void q50_density_blur(float& density, int4& q50_vox, const Ray& ray,
                                 float const energy_in, const int4& vox, const float max_step);

// Average blur
__device__ void get_average_blur_step(float& step,
                                      float& step_water,
                                      float& de,
                                      const float max_step,
                                      const Ray& ray,
                                      const int4& vox);
__device__ void average_blur (float& density, float& mass_stp_ratio,
                              const Ray& ray, const int4& vox, const float max_step);

__device__ void step_blur_get_voxels(int* voxels, const Ray& ray, const float& energy_in,
                                     const int4& vox, const int& n, const float& scale,
                                     const float max_step);


__device__ float get_residual_range (float const energy);

__device__ float get_energy_loss (float const step_water, float const energy);

__device__ float massStpRatio(const float energy, const int4& vox);

__device__ float get_dose_at (float const energy, float const wepl_r,
                               float const wepl_d);

#endif