#ifndef __RAY_CUH__
#define __RAY_CUH__

class Ray
{
public:
    float3 pos;
    float3 dir;
    float wepl;
    float energy;

    __device__ Ray(float4 x_, float4 vx_);
    __device__ float3 position();
    __device__ float3 direction();
    __device__ void step(float step, float density);
    __device__ bool isAlive();

private:
    bool alive;
};


#endif
