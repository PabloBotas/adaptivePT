#include "gpu_physics.cuh"

#include "gpu_device_globals.cuh"

__device__ void get_water_step(float& step,
                               float& step_water,
                               float& de,
                               const float max_step,
                               const float energy_in,
                               const int4& vox)
{
    // Get density
    float const density  = tex3D(dens_tex, vox.z, vox.y, vox.x);
    // Get stp ratio
    float const mass_stp_ratio = massStpRatio(energy_in, vox);

    // Set steps
    step = max_step;
    step_water = mass_stp_ratio*density*max_step;
    de = get_energy_loss (step_water, energy_in);
    // Verify it's not too much
    if (de >= energy_in)
    {
        de = energy_in;
        step_water = get_residual_range (energy_in);
        step = step_water/(mass_stp_ratio*density);
    }

    printf("Step: %f %f %f %f %f %f %f\n", step_water,
                                           step, max_step,
                                           density, mass_stp_ratio,
                                           de, energy_in);
}

__device__ float get_residual_range (float const energy)
// find eloss based on water track length s and initial energy
// Based on Kawrakow, 1999
{
    // MP, stp_w_delta_e, stp_w_min_e, stp_w_tex and stp_w_b_coeff_tex are globals
    float index = (energy/2 - stp_w_min_e)/stp_w_delta_e + 0.5f;
    float stop_power = tex1D(stp_w_tex, index);
    float b = tex1D(stp_w_b_coeff_tex, index);

    float tau = energy/(2*MP);

    float a = 1 +
        ( 2*b*b*(1+tau)*(1+tau)*(2+tau) + b*(1+tau)*(tau*tau+3*tau-2) - 6*tau ) /
        ( 6*(1+tau)*(1+tau)*(2+tau) );

    float step_water = energy/stop_power * a;

    return step_water;
}

__device__ float get_energy_loss (float const step_water, float const energy)
// find eloss based on water track length s and initial energy
// Based on Kawrakow, 1999
{
    // MP, stp_w_delta_e, stp_w_min_e, stp_w_tex and stp_w_b_coeff_tex are globals
    float index = (energy - stp_w_min_e)/stp_w_delta_e + 0.5f;
    float de1 = step_water*tex1D(stp_w_tex, index);
    float b = tex1D(stp_w_b_coeff_tex, index);

    float temp = energy/MP;
    float eps = de1/energy;

    float de = de1*( 1 +
               eps/(1+temp)/(2+temp) +
               eps*eps*(2+2*temp+temp*temp)/(1+temp)/(1+temp)/(2+temp)/(2+temp) -
               b*eps*(0.5+2*eps/3/(1+temp)/(2+temp) +
               (1-b)*eps/6) );

    de = fmin(de, energy);
    return de;
}

__device__ float massStpRatio(const float energy, const int4& vox)
//  mass stopping power ratio wrt water for a material
{
    float const energy_index = (energy - stp_ratio_min_e)/stp_ratio_delta_e + 0.5;
    int const material_id = tex3D(matid_tex, vox.z, vox.y, vox.x);
    int const material_index = material_id + 0.5;
    return tex2D(stp_ratio_tex, energy_index , material_index);
}