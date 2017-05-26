#include "gpu_ray_class.cuh"

__device__ Ray::Ray(float4 x_, float4 vx_, short2 ix_) : _alive(true)
{
    pos.x        = x_.x;
    pos.y        = x_.y;
    pos.z        = x_.z;
    initial_wepl = x_.w;
    dir.x        = vx_.x;
    dir.y        = vx_.y;
    dir.z        = vx_.z;
    energy       = vx_.w;

    wepl         = initial_wepl;

    beam_id      = ix_.x;
    spot_id      = ix_.y;
}

__device__ float3 Ray::position()
{
    return pos;
}

__device__ float3 Ray::direction()
{
    return dir;
}

__device__ void Ray::move(float step, float step_water)
{
    pos.x += step*dir.x;
    pos.y += step*dir.y;
    pos.z += step*dir.z;

    wepl -= step_water;
    if (wepl <= _min_wepl)
        _alive = false;
}

__device__ bool Ray::isAlive()
{
    return _alive;
}
