#include "gpu_ray_class.cuh"

__device__ Ray::Ray(float4 x_, float4 vx_) : alive(true)
{
    pos.x  = x_.x;
    pos.y  = x_.y;
    pos.z  = x_.z;
    wepl   = x_.w;
    dir.x  = vx_.x;
    dir.y  = vx_.y;
    dir.z  = vx_.z;
    energy = vx_.w;
}

__device__ float3 Ray::position()
{
    return pos;
}

__device__ float3 Ray::direction()
{
    return dir;
}

__device__ void Ray::step(float step, float density)
{
    float step_wepl = density*step;
    if(step_wepl >= wepl)
    {
        alive = false;
        step  = wepl/density;
    }
    pos.x += step*dir.x;
    pos.y += step*dir.y;
    pos.z += step*dir.z;
    wepl -= step_wepl;
}

__device__ bool Ray::isAlive()
{
    return alive;
}
