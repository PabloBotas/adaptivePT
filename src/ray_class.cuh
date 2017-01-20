#ifndef __RAY_CUH__
#define __RAY_CUH__

class Ray
{
public:
    float3 pos;
    float3 dir;

    __device__ Ray(float3 pos, float3 dir);
    __device__ float3 position();
    __device__ float3 direction();
    __device__ void step(float step);
};


#endif
