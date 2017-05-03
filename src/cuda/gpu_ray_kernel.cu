#include "gpu_ray_kernel.cuh"
#include "gpu_device_interaction.cuh"
#include "gpu_ray_class.cuh"
#include "gpu_geometry_operations.cuh"

__global__ void calculateRays_kernel(const int num,
		                             float4* endpoints,
		                             float* traces,
		                             const short* spots_per_beam)
{
    const int id = blockIdx.x*blockDim.x + threadIdx.x;
    if(id < num)
    {
        Ray ray(xdata[id], vxdata[id], ixdata[id]);
        int4 vox = make_int4(ray.pos/ctVoxSize);
        vox.w = getVoxelIndex(vox);

        VoxelUpdater voxUpdater;
        VoxelStepper voxStepper;
        if (id == 0)
        {
//        	printf("%d %d %d -> %f %f %f -> %f MeV -> %f -> %d %d\n", vox.x, vox.y, vox.z, ray.pos.x, ray.pos.y, ray.pos.z, ray.energy/1000000, ray.initial_wepl, ray.isAlive(), vox.w);
        	printf("ID - NX NY NZ -> X Y Z -> WEPL\n");
        }
        if (id < 10)
            printf("%d - %d %d %d -> %f %f %f -> %f\n", id, vox.x, vox.y, vox.z, ray.pos.x, ray.pos.y, ray.pos.z, ray.wepl);

        while (ray.isAlive())
        {
        	if (vox.w != -1)
        		atomicAdd(&traces[vox.w], 1.f);
        	float density = tex3D(dens_tex, vox.z, vox.y, vox.x);
            float step    = inters(ray, vox, voxUpdater, voxStepper);
            if (step == INF)
            	break;
//        	if (id == 0)
//        		printf("%d %d %d -> %f %f %f -> %f -> %f -> %f\n", vox.x, vox.y, vox.z, ray.pos.x, ray.pos.y, ray.pos.z, ray.wepl, step, density);
            ray.step(step, density);
            changeVoxel(vox, voxUpdater, voxStepper);
        }

        if (id < 10)
        	printf("%d - %d %d %d -> %f %f %f -> %f\n", id, vox.x, vox.y, vox.z, ray.pos.x, ray.pos.y, ray.pos.z, ray.wepl);
		atomicAdd(&traces[vox.w], 50);

        unsigned int index = get_endpoints_index(ray.beam_id, ray.spot_id, spots_per_beam);
        endpoints[index].x = ray.pos.x;
        endpoints[index].y = ray.pos.y;
        endpoints[index].z = ray.pos.z;
        endpoints[index].w = ray.wepl;
    }
}

__device__ unsigned int get_endpoints_index(const short beam_id,
                                            const short spot_id,
                                            const short* spots_per_beam)
{
    unsigned int index = spot_id;
    for (int ibeam = 0; ibeam < beam_id; ibeam++)
    {
        index += spots_per_beam[ibeam];
    }
    return index;
}

__global__ void rays_to_device_kernel(const int num,
                                      float2* angles,
                                      const float3 ct_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num)
    {
        float4 pos  = xdata[tid];
        float4 vel  = vxdata[tid];
        short2 meta = ixdata[tid]; // x = beam_id, y = spot_id

        //  get angles
        float phi   =   angles[meta.x].x; // gantry angle
        float theta = - angles[meta.x].y; // couch angle

        //  offset locations, the coordinate of corner A in internal coordinate

        float3 offset = ct_offsets - 0.5*ctVox*ctVoxSize;

        //  rotate location using gantry and couch
        pos.x = pos.x*__cosf(theta) - __sinf(theta)*(pos.y*__sinf(phi) + pos.z*__cosf(phi)) - offset.x;
        pos.y = (pos.y*__cosf(phi) - pos.z*__sinf(phi)) - offset.y;
        pos.z = pos.x*__sinf(theta) + __cosf(theta)*(pos.y*__sinf(phi) + pos.z*__cosf(phi)) - offset.z;
        xdata[tid] = pos;

        //  velocity
        vel.x = vel.x*__cosf(theta) - __sinf(theta)*(vel.y*__sinf(phi) + vel.z*__cosf(phi));
        vel.y = (vel.y*__cosf(phi) - vel.z*__sinf(phi));
        vel.z = vel.x*__sinf(theta) + __cosf(theta)*(vel.y*__sinf(phi) + vel.z*__cosf(phi));
        vxdata[tid] = vel;

        float temp;
        float alphaMin = -1.0f;
        float alpha1x, alphanx, alpha1y, alphany, alpha1z, alphanz;

        alpha1x = (0.1f* ctVoxSize.x - pos.x)/vel.x;
        alphanx = (ctVox.x * ctVoxSize.x - 0.1f* ctVoxSize.x - pos.x)/vel.x;
        alpha1y = (0.1f* ctVoxSize.y - pos.y)/vel.y;
        alphany = (ctVox.y * ctVoxSize.y - 0.1f* ctVoxSize.y - pos.y)/vel.y;
        alpha1z = (0.1f* ctVoxSize.z - pos.z)/vel.z;
        alphanz = (ctVox.z * ctVoxSize.z - 0.1f* ctVoxSize.z - pos.z)/vel.z;

        if((alpha1x<0.0f && alphanx<0.0f) ||(alpha1y<0.0f && alphany<0.0f)||(alpha1z<0.0f && alphanz<0.0f))
        {

        }
        else if((alpha1x*alphanx<=0.0f) && (alpha1y*alphany<=0.0f) && (alpha1z*alphanz<=0.0f))
        {

        }
        else
        {
            temp = min(alpha1x, alphanx);
            alphaMin = max(alphaMin, temp);

            temp = min(alpha1y, alphany);
            alphaMin = max(alphaMin, temp);

            temp = min(alpha1z, alphanz);
            alphaMin = max(alphaMin, temp);

            xdata[tid].x = pos.x + vel.x*alphaMin;
            xdata[tid].y = pos.y + vel.y*alphaMin;
            xdata[tid].z = pos.z + vel.z*alphaMin;
        }
    }
}
