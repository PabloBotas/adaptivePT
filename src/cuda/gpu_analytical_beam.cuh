#ifndef __ANALYTICAL_BEAM_CUH__
#define __ANALYTICAL_BEAM_CUH__


class AnalyticalBeam
{
    float b;
    float n;
    float s;
    float w;
public:
    // CONSTRUCTORS ----------------------------------------------------------
    __device__ AnalyticalBeam (float const energy, float const wepl_d);
    __device__ ~AnalyticalBeam ();
    __device__ float get_dose_at (float const wepl_r);
    __device__ void get_pars (float pars[]);
    __device__ float get_sigma ();
};


#endif
