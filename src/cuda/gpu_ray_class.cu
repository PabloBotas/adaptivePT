#include "gpu_ray_class.cuh"

#include "gpu_device_globals.cuh"
#include "gpu_physics.cuh"
#include "helper_math.h"

// CONSTRUCTORS --------------------------------------------------------------
__device__ Ray::Ray (double4 x_, double4 vx_, short2 ix_)
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

// BASIC GETTERS -------------------------------------------------------------

__device__ double3 Ray::get_position ()
{
    return pos;
}

__device__ double3 Ray::get_direction ()
{
    return dir;
}

__device__ double Ray::get_wepl ()
{
    return wepl;
}

__device__ double Ray::get_energy ()
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

// GETTERS -------------------------------------------------------------------

__device__ bool Ray::is_alive()
{
    _alive = energy > 0;
    return _alive;
}

// BASIC SETTERS -------------------------------------------------------------

__device__ void Ray::set_energy(double m)
{
    energy = m;
}

__device__ void Ray::set_wepl(double m)
{
    wepl = m;
}

__device__ void Ray::set_position (double4 d)
{
    set_position(make_double3(d));
}

__device__ void Ray::set_position (double3 d)
{
    pos = d;
}

__device__ void Ray::set_direction (double4 d)
{
    set_direction(make_double3(d));
}

__device__ void Ray::set_direction (double3 d)
{
    dir = normalize(d);
}

// SETTERS -------------------------------------------------------------------

__device__ void Ray::set_direction_to_point (double4 p)
{
    double3 p2 = make_double3(p);
    set_direction_to_point(p2);
}

__device__ void Ray::set_direction_to_point (double3 p)
{
    dir = p - pos;
    double norm = length(dir);
    dir /= norm;
}

// ACTIONS -------------------------------------------------------------------

__device__ void Ray::move (const double& step,
                           const double& step_water,
                           const double& de)
{
    pos += step*dir;
    energy -= de;
    wepl += step_water;
    if (energy <= stp_w_min_e)
        _alive = false;
}

__device__ void Ray::kill()
{
    energy = 0;
    _alive = false;
}

__device__ void Ray::print()
{
    printf("Beam ID:   %d\n"
           "Spot ID:   %d\n"
           "Position:  %f %f %f\n"
           "Direction: %f %f %f\n"
           "Energy:    %f\n"
           "WEPL:      %f\n",
           beam_id, spot_id, pos.x, pos.y, pos.z,
           dir.x, dir.y, dir.z, energy, wepl);
}

