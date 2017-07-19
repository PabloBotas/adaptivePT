#ifndef __RAY_CUH__
#define __RAY_CUH__

class Ray
{
    double wepl;
    double energy;
    short  beam_id;
    short  spot_id;
public:
    double3 pos;
    double3 dir;

    __device__ Ray (double4 x_, double4 vx_, short2 ix);

//    __device__ double energyToRange (double energy);
    __device__ double3 get_position ();
    __device__ double3 get_direction ();
    __device__ short get_beam_id ();
    __device__ short get_spot_id ();
    __device__ double get_energy ();
    __device__ double get_wepl ();
    __device__ void move (const double& step, const double& wepl, const double& de);
    __device__ bool is_alive ();
    __device__ void set_energy (double e);
    __device__ void set_wepl (double e);
    __device__ void set_direction_to_point (double3 d);
    __device__ void set_direction_to_point (double4 d);
    __device__ void set_position (double4 d);
    __device__ void set_position (double3 d);
    __device__ void print();
    __device__ void kill();


private:
    bool  _alive    = true;
    double _min_wepl = 1e-08;
};


#endif
