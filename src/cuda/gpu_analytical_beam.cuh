#ifndef __ANALYTICAL_BEAM_CUH__
#define __ANALYTICAL_BEAM_CUH__


class AnalyticalBeam
{
    double b;
    double n;
    double s;
    double w;
public:
    // CONSTRUCTORS ----------------------------------------------------------
    __device__ AnalyticalBeam (double const energy, double const wepl_d);
    __device__ ~AnalyticalBeam ();
    __device__ double get_dose_at (double const wepl_r);
    __device__ void get_pars (double pars[]);
    __device__ double get_sigma ();
};


#endif
