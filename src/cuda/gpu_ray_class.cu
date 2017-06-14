#include "gpu_ray_class.cuh"

#include "gpu_device_globals.cuh"
#include "helper_math.h"

__device__ Ray::Ray (float4 x_, float4 vx_, short2 ix_)
{
    pos.x  = x_.x;
    pos.y  = x_.y;
    pos.z  = x_.z;
    wepl   = x_.w;
    dir.x  = vx_.x;
    dir.y  = vx_.y;
    dir.z  = vx_.z;
    energy = vx_.w;

    beam_id = ix_.x;
    spot_id = ix_.y;
}

__device__ float3 Ray::get_position ()
{
    return pos;
}

__device__ float3 Ray::get_direction ()
{
    return dir;
}

__device__ float Ray::get_wepl ()
{
    return wepl;
}

__device__ float Ray::get_energy ()
{
    return energy;
}

__device__ short Ray::get_beam_id ()
{
    return beam_id;
}

__device__ short Ray::get_spot_id ()
{
    return spot_id;
}

__device__ void Ray::reverse_direction ()
{
    dir.x = -dir.x;
    dir.y = -dir.y;
    dir.z = -dir.z;
}

__device__ void Ray::set_direction (float4 p)
{
    float3 p2 = make_float3(p);
    set_direction(p2);
}

__device__ void Ray::set_direction (float3 p)
{
    dir = p - pos;
    float norm = length(dir);
    dir /= norm;
}

__device__ void Ray::move (const float step, const float step_water)
{
    pos += step*dir;

    lose_energy(step_water);
    wepl -= step_water;
    if (energy <= stp_w_min_e || wepl <= _min_wepl)
        _alive = false;
}

__device__ void Ray::lose_energy (float const step_water)
// find eloss based on water track length s and initial energy
// Based on Kawrakow, 1999
{
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

    de = (de < energy) ? de : energy;
    energy -= de;
}

__device__ void Ray::set_energy(float m)
{
    energy = m;
}

__device__ bool Ray::is_alive()
{
    return _alive;
}
