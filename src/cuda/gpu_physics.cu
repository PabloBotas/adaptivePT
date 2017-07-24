#include "gpu_physics.cuh"

#include "gpu_device_globals.cuh"

__device__ void get_step(double& step,
                         double& step_water,
                         double& de,
                         const double max_step,
                         const double energy_in,
                         const int4& vox)
{
    // Get density
    double const density  = tex3D(dens_tex, vox.z, vox.y, vox.x);
    // Get stp ratio
    double const mass_stp_ratio = massStpRatio(energy_in, vox);
    // Set steps
    step = max_step;
    step_water = mass_stp_ratio*density*max_step;
    de = get_energy_loss (step_water, energy_in);

    if (de == energy_in)
    {
        step_water = get_residual_range (energy_in);
        step = step_water/(massStpRatio(energy_in/2, vox)*density);
    }
}

__device__ double get_residual_range (double const energy)
{
    // Find range based on half the residual energy,
    // the get some kind of average stopping power

    // Function based on Kawrakow (1999):
    // Accurate condense history Monte Carlo simulation of electron transport.
    // I. EGSnrc, the new EGS4 version

    // MP, stp_w_delta_e, stp_w_min_e, stp_w_tex and stp_w_b_coeff_tex are globals
    float index = (energy/2 - stp_w_min_e)/stp_w_delta_e + 0.5f;
    double stop_power = tex1D(stp_w_tex, index);
    double b = tex1D(stp_w_b_coeff_tex, index);

    double tau = energy/(2*MP);

    double a = 1 +
             ( 2*b*b*(1+tau)*(1+tau)*(2+tau) + b*(1+tau)*(tau*tau+3*tau-2) - 6*tau ) /
             ( 6*(1+tau)*(1+tau)*(2+tau) );

    double step_water = energy/stop_power * a;

    return step_water;
}

__device__ double get_energy_loss (double const step_water, double const energy)
// find eloss based on water track length s and initial energy
// Based on Kawrakow, 1999
{
    // MP, stp_w_delta_e, stp_w_min_e, stp_w_tex and stp_w_b_coeff_tex are globals
    float index = (energy - stp_w_min_e)/stp_w_delta_e + 0.5f;
    double de1 = step_water*tex1D(stp_w_tex, index);
    double b = tex1D(stp_w_b_coeff_tex, index);

    double tau = energy/MP;
    double eps = de1/energy;

    double de = de1*( 1 +
                eps/(1+tau)/(2+tau) +
                eps*eps*(2+2*tau+tau*tau)/(1+tau)/(1+tau)/(2+tau)/(2+tau) -
                b*eps*(0.5+2*eps/3/(1+tau)/(2+tau) +
                (1-b)*eps/6) );

    de = fmin(de, energy);
    return de;
}

__device__ double massStpRatio(const double energy, const int4& vox)
//  mass stopping power ratio wrt water for a material
{
    float const energy_index = (energy - stp_ratio_min_e)/stp_ratio_delta_e + 0.5;
    int const material_id = tex3D(matid_tex, vox.z, vox.y, vox.x);
    int const material_index = material_id + 0.5;
    return tex2D(stp_ratio_tex, energy_index , material_index);
}
