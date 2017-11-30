#include "gpu_analytical_beam.cuh"

#include "gpu_device_globals.cuh"
#include "helper_math.h"


// CONSTRUCTORS --------------------------------------------------------------
__device__ AnalyticalBeam::AnalyticalBeam (float const energy, float const wepl_d)
{
    float energy_index = (energy - bp_energy_min)/bp_energy_delta + 0.5;
    float range_fraction = wepl_d/tex1D(bp_range_tex, energy_index);

    float energy_index_0 = floor(energy_index - 0.5) + 0.5;
    float energy_index_1 = energy_index_0 + 1;
    float equiv_wepl_0 = range_fraction*tex1D(bp_range_tex, energy_index_0);
    float equiv_wepl_1 = range_fraction*tex1D(bp_range_tex, energy_index_1);
    float equiv_depth_index_0 = equiv_wepl_0/bp_depth_delta + 0.5;
    float equiv_depth_index_1 = equiv_wepl_1/bp_depth_delta + 0.5;

    float b_0 = tex2D(bp_b_tex, energy_index_0, equiv_depth_index_0);
    float n_0 = tex2D(bp_n_tex, energy_index_0, equiv_depth_index_0);
    float s_0 = tex2D(bp_s_tex, energy_index_0, equiv_depth_index_0);
    float w_0 = tex2D(bp_w_tex, energy_index_0, equiv_depth_index_0);
    float b_1 = tex2D(bp_b_tex, energy_index_1, equiv_depth_index_1);
    float n_1 = tex2D(bp_n_tex, energy_index_1, equiv_depth_index_1);
    float s_1 = tex2D(bp_s_tex, energy_index_1, equiv_depth_index_1);
    float w_1 = tex2D(bp_w_tex, energy_index_1, equiv_depth_index_1);

    float weight = energy_index-energy_index_0;
    b = (1-weight)*b_0 + weight*b_1;
    n = (1-weight)*n_0 + weight*n_1;
    s = (1-weight)*s_0 + weight*s_1;
    w = (1-weight)*w_0 + weight*w_1;
}


__device__ AnalyticalBeam::~AnalyticalBeam ()
{
}

// BASIC GETTERS -------------------------------------------------------------
__device__ float AnalyticalBeam::get_dose_at (float const R)
{
    if (n == 0) {
        return 0;
    } else {
        float norma = w/(sqrtf(2*PI)*s) + 2*(1-w)/(PI*sqrtf(b));
        float gauss = 1/(sqrtf(2*PI)*s) * exp(-0.5 * R*R/(s*s));
        float ruthe = 2*pow(b, 1.5)/(PI*pow(R*R + b, 2));
        return fabs(n/norma * (w*gauss + (1-w)*ruthe));
    }
    
}


__device__ void AnalyticalBeam::get_pars (float pars[])
{
    pars[0] = n;
    pars[1] = w;
    pars[2] = s;
    pars[3] = b;
}


__device__ float AnalyticalBeam::get_sigma ()
{
    return s;
}
