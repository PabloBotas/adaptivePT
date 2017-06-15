#ifndef __RAY_CUH__
#define __RAY_CUH__

class Ray
{
    float wepl;
    float energy;
    short beam_id;
    short spot_id;
public:
    float3 pos;
    float3 dir;

    __device__ Ray (float4 x_, float4 vx_, short2 ix);

//    __device__ float energyToRange (float energy);
    __device__ float3 get_position ();
    __device__ float3 get_direction ();
    __device__ short get_beam_id ();
    __device__ short get_spot_id ();
    __device__ float get_energy ();
    __device__ float get_wepl ();
    __device__ void lose_energy (float const step);
    __device__ void reverse_direction ();
    __device__ void move (const float step, const float wepl);
    __device__ bool is_alive ();
    __device__ void set_energy (float e);
    __device__ void set_wepl (float e);
    __device__ void set_direction (float3 d);
    __device__ void set_direction (float4 d);
    __device__ void set_position (float4 d);
    __device__ void set_position (float3 d);
    __device__ void print();
    __device__ void kill();


private:
    bool  _alive    = true;
    float _min_wepl = 1e-08;
};


#endif
