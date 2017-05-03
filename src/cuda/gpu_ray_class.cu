#include "gpu_ray_class.cuh"

__device__ Ray::Ray(float4 x_, float4 vx_, short2 ix_) : alive(true)
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

__device__ void Ray::step(float step, float density)
{
    float step_wepl = density*step;
    if(wepl < step_wepl)
	{
    	step_wepl = wepl;
		alive = false;
		step = step_wepl/density;
	}
    wepl -= step_wepl;

    pos.x += step*dir.x;
    pos.y += step*dir.y;
    pos.z += step*dir.z;
}

__device__ bool Ray::isAlive()
{
    return alive;
}

//__device__ float Ray::energyToRange(float energy)
//{
//	// Foundation of an analytical proton beamlet model for inclusion
//	// in a general proton dose calculation system (sup. materials)
//    // Ulmer et al., Rad. Phys. and Chem. 80(2011)378
//	const short N = 4;
//	float a[N] = {6.94656e-3, 8.13116e-4, -1.21068e-6, 1.053e-9};
//
//    float range = 0;
//    energy *= 1e-6; // convert to MeV
//    for(int i = 1; i <= N; i++)
//    	range += a[i-1]*powf(energy, i);
//
//    return range;
//}
