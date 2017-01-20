
#include "ray_class.cuh"

__device__ Ray::Ray(float3 pos_, float3 dir_)
{
    pos.x = pos_.x;
    pos.y = pos_.y;
    pos.z = pos_.z;
    dir.x = dir_.x;
    dir.y = dir_.y;
    dir.z = dir_.z;
}

__device__ float3 Ray::position()
{
    return pos;
}

__device__ float3 Ray::direction()
{
    return dir;
}

__device__ void Ray::step(float step)
{
   pos.x += step*dir.x;
   pos.y += step*dir.y;
   pos.z += step*dir.z;
}
