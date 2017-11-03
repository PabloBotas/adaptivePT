#include "gpu_analytical_beam.cuh"

#include "gpu_device_globals.cuh"
#include "helper_math.h"


// CONSTRUCTORS --------------------------------------------------------------
__device__ AnalyticalBeam::AnalyticalBeam (double const energy, double const wepl_d)
{
    energy_index = (energy - bp_energy_min)/bp_energy_delta + 0.5;
    depth_index = wepl_d/bp_depth_delta + 0.5;
    fill_pars();
}


__device__ AnalyticalBeam::~AnalyticalBeam ()
{
}

// BASIC GETTERS -------------------------------------------------------------


__device__ void AnalyticalBeam::fill_pars ()
{
    n = tex2D(bp_n_tex, energy_index, depth_index);
    if (n > 0) {
        b = tex2D(bp_b_tex, energy_index, depth_index);
        s = tex2D(bp_s_tex, energy_index, depth_index);
        w = tex2D(bp_w_tex, energy_index, depth_index);
    } else {
        b = 1;
        s = 1;
        w = 1;
    }
    
}


__device__ double AnalyticalBeam::get_dose_at (double const wepl_r)
{
    if (n == 0)
        return 0;
    else
        return fabs(n * (w/(sqrtf(2*PI)*s) * exp(-0.5 * wepl_r*wepl_r/(s*s)) +
               (1-w) * 2*pow(b, 1.5)/PI *1/pow(wepl_r*wepl_r + b, 2)) );
}


__device__ void AnalyticalBeam::get_pars (double pars[])
{
    pars[0] = n;
    pars[1] = w;
    pars[2] = s;
    pars[3] = b;
}


__device__ double AnalyticalBeam::get_sigma ()
{
    return s;
}
