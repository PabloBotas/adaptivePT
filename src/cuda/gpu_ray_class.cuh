#ifndef __RAY_CUH__
#define __RAY_CUH__

class Ray
{
public:
    float3 pos;
    float3 dir;
    float wepl;
    float initial_wepl;
    float energy;
    short beam_id;
    short spot_id;

    __device__ Ray(float4 x_, float4 vx_, short2 ix);
//    __device__ float energyToRange(float energy);
    __device__ float3 position();
    __device__ float3 direction();
    __device__ void move(float step, float wepl);
    __device__ bool isAlive();

private:
    bool _alive;
    float _min_wepl = 1e-08;
};


#endif
