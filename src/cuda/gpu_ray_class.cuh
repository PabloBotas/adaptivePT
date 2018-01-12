#ifndef __RAY_CUH__
#define __RAY_CUH__

class Ray
{
    float wepl;
    float energy;
    float initial_energy;
    short beam_id;
    short spot_id;
public:
    float3 pos;
    float3 dir;

    // CONSTRUCTORS ----------------------------------------------------------
    __device__ Ray (float4 x_, float4 vx_, short2 ix);

    // BASIC GETTERS ---------------------------------------------------------
    __device__ float3 get_position () const;
    __device__ float3 get_direction () const;
    __device__ short get_beam_id () const;
    __device__ short get_spot_id () const;
    __device__ float get_energy () const;
    __device__ float get_initial_energy () const;
    __device__ float get_wepl () const;
    // GETTERS ---------------------------------------------------------------
    __device__ bool is_alive ();
    // BASIC SETTERS ---------------------------------------------------------
    __device__ void set_energy (float e);
    __device__ void set_wepl (float e);
    __device__ void set_position (float4 d);
    __device__ void set_position (float3 d);
    __device__ void set_direction (float4 d);
    __device__ void set_direction (float3 d);
    // SETTERS ---------------------------------------------------------------
    __device__ void set_direction_to_point (float3 d);
    __device__ void set_direction_to_point (float4 d);
    // ACTIONS ---------------------------------------------------------------    
    __device__ void move (const float& step, const float& wepl, const float& de);
    __device__ void kill();
    // QUERIES ---------------------------------------------------------------
    __device__ void print();


private:
    bool  _alive     = true;
    float _min_wepl = 1e-08;
    float error_wepl = 0.f;
    float error_energy = 0.f;
    float3 error_pos = make_float3(0.f, 0.f, 0.f);
    float3 error_dir = make_float3(0.f, 0.f, 0.f);
};

#endif
