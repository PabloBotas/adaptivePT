#include "gpu_physics.cuh"

#include "gpu_device_globals.cuh"
#include "gpu_ray_class.cuh"
#include "gpu_analytical_beam.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_utils.cuh"

#include <curand.h>
#include <curand_kernel.h>


__device__ void get_step(float& step,
                         float& step_water,
                         float& de,
                         const float max_step,
                         const Ray& ray,
                         const int4& vox)
{
#if defined(__STEP_CENTRAL_AXIS__)
    get_step(step, step_water, de, max_step, ray.get_energy(), vox);
#elif defined(__STEP_Q50__)
    get_q50_blur_step(step, step_water, de, max_step, ray, vox);
#else
    get_average_blur_step(step, step_water, de, max_step, ray, vox);
#endif
}


__device__ void get_step(float& step,
                         float& step_water,
                         float& de,
                         const float max_step,
                         const float energy_in,
                         const int4& vox)
{
    // Get density
    float const density = tex3D(dens_tex, vox.z, vox.y, vox.x);
    // Get stp ratio
    float const mass_stp_ratio = massStpRatio(energy_in, vox);
    // Set steps
    step = max_step;
    step_water = mass_stp_ratio*density*max_step;
    de = get_energy_loss (step_water, energy_in);
    // const int thread = blockIdx.x*blockDim.x + threadIdx.x;
    // if (thread == 0) {
    //     printf("%f %f ", density, mass_stp_ratio);
    // }

    if (de > energy_in) {
        step_water = get_residual_range (energy_in);
        step = step_water/(massStpRatio(energy_in/2, vox)*density);
    }
}


__device__ void get_q50_blur_step(float& step,
                                  float& step_water,
                                  float& de,
                                  const float max_step,
                                  const Ray& ray,
                                  const int4& vox)
{
    // Select density and voxel corresponding to the Q50 across a 9 position
    float ray_energy = ray.get_energy();
    int4 q50_vox;
    float density=0;
    q50_density_blur(density, q50_vox, ray, ray.get_initial_energy(), vox, max_step);

    // Get stp ratio
    float const mass_stp_ratio = massStpRatio(ray_energy, q50_vox);
    
    // Set steps
    step = max_step;
    step_water = mass_stp_ratio*density*max_step;
    de = get_energy_loss (step_water, ray_energy);

    if (de > ray_energy) {
        step_water = get_residual_range (ray_energy);
        step = step_water/(massStpRatio(ray_energy/2, q50_vox)*density);
    }
}


__device__ void q50_density_blur(float& density, int4& q50_vox, const Ray& ray,
                                 float const energy_in, const int4& vox, const float max_step)
{
    int const N = 9;
    int const q50_index = int(N/2+0.5);
    float const sigma_scale = 0.5;
    int voxels_array[N];
    // Set 8 points at a distance from the central axis
    step_blur_get_voxels(voxels_array, ray, energy_in, vox, N, sigma_scale, max_step);
    // Calculate 50 percentile of the density
    float density_array[N];
    for (int i = 0; i < N; ++i) {
        int4 voxid = make_int4(getVoxelCoords(voxels_array[i]), voxels_array[i]);
        density_array[i] = tex3D(dens_tex, voxid.z, voxid.y, voxid.x);
    }
    bubble_sort(density_array, N, voxels_array);
    density = density_array[q50_index]; // Approximately the 50% percentile
    // Select the voxel showing that density
    q50_vox = make_int4(getVoxelCoords(voxels_array[q50_index]), voxels_array[q50_index]);
}


__device__ void get_average_blur_step(float& step,
                                      float& step_water,
                                      float& de,
                                      const float max_step,
                                      const Ray& ray,
                                      const int4& vox)
{
    // Select density and voxel corresponding to the Q50 across a 9 position
    float density=0, mass_stp_ratio=0;
    average_blur(density, mass_stp_ratio, ray, vox, max_step);
    float ray_energy = ray.get_energy();
    // Set steps
    step = max_step;
    step_water = mass_stp_ratio*density*max_step;
    de = get_energy_loss (step_water, ray_energy);

    if (de > ray_energy) {
        step_water = get_residual_range (ray_energy);
        step = step_water/(massStpRatio(ray_energy/2, vox)*density);
    }
}


__device__ void average_blur (float& density, float& mass_stp_ratio,
                              const Ray& ray, const int4& vox, const float max_step)
{
    int const N = 9;
    float sigma_scale = 0; // this value is never used
#if defined(__STEP_AVERAGE__)
    sigma_scale = 1;
#elif defined(__STEP_AVERAGE2__)
    sigma_scale = 0.5;
#elif defined(__STEP_AVERAGE3__)
    sigma_scale = 0.25;
#endif
    int voxels_array[N];
    // Set 8 points at a distance from the central axis
    step_blur_get_voxels(voxels_array, ray, ray.get_initial_energy(), vox, N, sigma_scale, max_step);
    float density_array[N], mass_stp_ratio_array[N];
    for (int i = 0; i < N; ++i) {
        int4 tempvox = make_int4(getVoxelCoords(voxels_array[i]), voxels_array[i]);
        density_array[i] = tex3D(dens_tex, tempvox.z, tempvox.y, tempvox.x);
        mass_stp_ratio_array[i] = massStpRatio(ray.get_energy(), tempvox);
    }

    // Find resultant density and STP
    float rel_weight = exp(-0.5*sigma_scale*sigma_scale);
    for (int i = 0; i < N-1; ++i) {
        density += rel_weight*density_array[i];
        mass_stp_ratio += rel_weight*mass_stp_ratio_array[i];
    }

    density += density_array[N-1];
    mass_stp_ratio += mass_stp_ratio_array[N-1];
    density /= (rel_weight*(N-1)+1);
    mass_stp_ratio /= (rel_weight*(N-1)+1);
}


__device__ void step_blur_get_voxels(int* voxels, const Ray& ray, const float& energy_in,
                                     const int4& vox, const int& n, const float& scale,
                                     const float max_step)
{
    // CALCULATE AND APPLY BLURRING AT BEGINNING OF STEP
    AnalyticalBeam beam(energy_in, ray.get_wepl());
    float sigma = scale*beam.get_sigma();

    // Get perpendicular vectors on ray plane v1, v2
    float3 u = ray.get_direction();
    float3 p = ray.get_position()+max_step/2;
    float3 v1, v2;
    if (u.z != 0)
        v1 = make_float3(1, 1, -(u.x+u.y)/u.z);
    else if (u.y != 0)
        v1 = make_float3(1, -(u.x+u.z)/u.y, 1);
    else
        v1 = make_float3(-(u.y+u.z)/u.x, 1, 1);
    v1 = normalize(v1);

    v2 = make_float3(u.y*v1.z-u.z*v1.y, u.x*v1.z-u.z*v1.x, u.x*v1.y-u.y*v1.x);
    v2 = normalize(v2);

    // Get a random number depending on the 3D position (unique)
    // The length of p is cast to an integer and close positions may
    // have the same integer values, so I multiply by 10000 to separate them
    curandState localState;
    curand_init(10000*length(p), 0, 0, &localState);
    float rand = 2*PI*curand_uniform(&localState);

    float3 temp;
    for (int i = 0; i < n-1; ++i) {
        // get positions around nominal ray and voxid
        temp = p + sigma*sin(i*2*PI/(n-1)+rand)*v1 + sigma*cos(i*2*PI/(n-1)+rand)*v2;
        float tempabs = get_voxel_abs(temp);
        if (tempabs < 0)
            tempabs = vox.w;
        voxels[i] = tempabs;
    }
    voxels[n-1] = vox.w;
}



__device__ float get_residual_range (float const energy)
{
    // Find range based on half the residual energy,
    // the get some kind of average stopping power

    // Function based on Kawrakow (1999):
    // Accurate condense history Monte Carlo simulation of electron transport.
    // I. EGSnrc, the new EGS4 version

    // MP, stp_w_delta_e, stp_w_min_e, stp_w_tex and stp_w_b_coeff_tex are globals
    float index = (energy/2 - stp_w_min_e)/stp_w_delta_e + 0.5f;
    float stop_power = tex1D(stp_w_tex, index);
    float tau = energy/(2*MP);
    float b = tex1D(stp_w_b_coeff_tex, index);
    float a = 1 +
             ( 2*b*b*(1+tau)*(1+tau)*(2+tau) + b*(1+tau)*(tau*tau+3*tau-2) - 6*tau ) /
             ( 6*(1+tau)*(1+tau)*(2+tau) );
    float step_water = energy/stop_power * a;

    return step_water;
}


__device__ float get_dose_at (float const energy, float const wepl_d, float const wepl_r)
{
    if (wepl_d < 0 || wepl_r < 0)
        return 0;
    AnalyticalBeam d(energy, wepl_d);
    return d.get_dose_at(wepl_r);
}


__device__ float get_energy_loss (float const step_water, float const energy)
// find eloss based on water track length s and initial energy
// Based on Kawrakow, 1999
{
    // MP, stp_w_delta_e, stp_w_min_e, stp_w_tex and stp_w_b_coeff_tex are globals
    float index = (energy - stp_w_min_e)/stp_w_delta_e + 0.5f;
    float de1 = step_water*tex1D(stp_w_tex, index);
    float b = tex1D(stp_w_b_coeff_tex, index);
    float tau = energy/MP;
    float eps = de1/energy;

    float de = de1*( 1 +
                eps/(1+tau)/(2+tau) +
                eps*eps*(2+2*tau+tau*tau)/(1+tau)/(1+tau)/(2+tau)/(2+tau) -
                b*eps*(0.5+2*eps/3/(1+tau)/(2+tau) +
                (1-b)*eps/6) );

    de = fmin(de, energy);
    return de;
}

__device__ float massStpRatio(const float energy, const int4& vox)
//  mass stopping power ratio wrt water for a material
{
    int const material_id = tex3D(matid_tex, vox.z, vox.y, vox.x);
    float const energy_index = (energy - stp_ratio_min_e)/stp_ratio_delta_e + 0.5;
    return tex2D(stp_ratio_tex, energy_index , material_id + 0.5);
}
