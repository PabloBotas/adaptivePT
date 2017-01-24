#include "gpu_ray_kernel.cuh"
#include "gpu_device_interaction.cuh"
#include "gpu_ray_class.cuh"
#include "gpu_geometry_operations.cuh"

__global__ void calculateRays_kernel(int num, float3 *scorer)
{
    const int id = blockIdx.x*blockDim.x + threadIdx.x;
    if(id < num)
    {
        Ray ray(xdata[id], vxdata[id]);
        int4 vox = make_int4(ray.pos/ctVoxSize);
        vox.w = getVoxelIndex(vox);

        VoxelUpdater voxUpdater;
        VoxelStepper voxStepper;
        while (ray.isAlive())
        {
            float density = tex3D(dens_tex, vox.z, vox.y, vox.x);
            float step    = inters(ray, vox, voxUpdater, voxStepper);
            ray.step(step, density);
            changeVoxel(vox, voxUpdater, voxStepper);
        }
        scorer[id] = ray.position();
    }
}


__global__ void rays_to_device_kernel(int num, float2 angles, float3 ct_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num)
    {
        //  get angles
        float phi   =   angles.x; // gantry angle
        float theta = - angles.y; // couch angle

        //  offset locations, the coordinate of corner A in internal coordinate

        float3 offset = ct_offsets - 0.5*ctVox*ctVoxSize;

        //  rotate location using gantry and couch
        float3 temp = make_float3(xdata[tid].x, xdata[tid].y, xdata[tid].z);
        float xs = temp.x*__cosf(theta) - __sinf(theta)*(temp.y*__sinf(phi) + temp.z*__cosf(phi)) - offset.x;
        float ys = (temp.y*__cosf(phi) - temp.z*__sinf(phi)) - offset.y;
        float zs = temp.x*__sinf(theta) + __cosf(theta)*(temp.y*__sinf(phi) + temp.z*__cosf(phi)) - offset.z;
        xdata[tid] = make_float4(xs, ys, zs, xdata[tid].w);

        //  velocity
        temp = make_float3(vxdata[tid].x, vxdata[tid].y, vxdata[tid].z);
        float vxs = temp.x*__cosf(theta) - __sinf(theta)*(temp.y*__sinf(phi) + temp.z*__cosf(phi));
        float vys = (temp.y*__cosf(phi) - temp.z*__sinf(phi));
        float vzs = temp.x*__sinf(theta) + __cosf(theta)*(temp.y*__sinf(phi) + temp.z*__cosf(phi));
        vxdata[tid] = make_float4(vxs, vys, vzs, vxdata[tid].w);

        float temp2;
        float alphaMin = -1.0f;
        float alpha1x, alphanx, alpha1y, alphany, alpha1z, alphanz;

        alpha1x = (0.1f* ctVoxSize.x - xs)/vxs;
        alphanx = (ctVox.x * ctVoxSize.x - 0.1f* ctVoxSize.x - xs)/vxs;
        alpha1y = (0.1f* ctVoxSize.y - ys)/vys;
        alphany = (ctVox.y * ctVoxSize.y - 0.1f* ctVoxSize.y - ys)/vys;
        alpha1z = (0.1f* ctVoxSize.z - zs)/vzs;
        alphanz = (ctVox.z * ctVoxSize.z - 0.1f* ctVoxSize.z - zs)/vzs;

        if((alpha1x<0.0f && alphanx<0.0f) ||(alpha1y<0.0f && alphany<0.0f)||(alpha1z<0.0f && alphanz<0.0f))
        {

        }
        else if((alpha1x*alphanx<=0.0f) && (alpha1y*alphany<=0.0f) && (alpha1z*alphanz<=0.0f))
        {

        }
        else
        {
            temp2 = min(alpha1x, alphanx);
            alphaMin = max(alphaMin, temp2);

            temp2 = min(alpha1y, alphany);
            alphaMin = max(alphaMin, temp2);

            temp2 = min(alpha1z, alphanz);
            alphaMin = max(alphaMin, temp2);

            xdata[tid].x = xs + vxs * alphaMin;
            xdata[tid].y = ys + vys * alphaMin;
            xdata[tid].z = zs + vzs * alphaMin;
        }
    }
}
