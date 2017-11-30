#ifndef __RAY_CUH__
#define __RAY_CUH__

class Ray
{
    double wepl;
    double energy;
    double initial_energy;
    short  beam_id;
    short  spot_id;
public:
    double3 pos;
    double3 dir;

    // CONSTRUCTORS ----------------------------------------------------------
    __device__ Ray (double4 x_, double4 vx_, short2 ix);

    // BASIC GETTERS ---------------------------------------------------------
    __device__ double3 get_position () const;
    __device__ double3 get_direction () const;
    __device__ short get_beam_id () const;
    __device__ short get_spot_id () const;
    __device__ double get_energy () const;
    __device__ double get_initial_energy () const;
    __device__ double get_wepl () const;
    // GETTERS ---------------------------------------------------------------
    __device__ bool is_alive ();
    // BASIC SETTERS ---------------------------------------------------------
    __device__ void set_energy (double e);
    __device__ void set_wepl (double e);
    __device__ void set_position (double4 d);
    __device__ void set_position (double3 d);
    __device__ void set_direction (double4 d);
    __device__ void set_direction (double3 d);
    // SETTERS ---------------------------------------------------------------
    __device__ void set_direction_to_point (double3 d);
    __device__ void set_direction_to_point (double4 d);
    // ACTIONS ---------------------------------------------------------------    
    __device__ void move (const double& step, const double& wepl, const double& de);
    __device__ void kill();
    // QUERIES ---------------------------------------------------------------
    __device__ void print();


private:
    bool  _alive     = true;
    double _min_wepl = 1e-08;
};


#endif
