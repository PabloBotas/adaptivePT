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
        int4 vox;
        vox.x = floor(ray.pos.x/ctVoxSize.x);
        vox.y = floor(ray.pos.y/ctVoxSize.y);
        vox.z = floor(ray.pos.z/ctVoxSize.z);
        vox.w = getVoxelIndex(vox);

        VoxelUpdater voxUpdater;
        VoxelStepper voxStepper;

        while (ray.isAlive() && vox.w != -1)
        {
            atomicAdd(&traces[vox.w], 1.f);
            float density = tex3D(dens_tex, vox.z, vox.y, vox.x);
            float step    = inters(ray, vox, voxUpdater, voxStepper);
            ray.step(step, density);
            changeVoxel(vox, voxUpdater, voxStepper);
        }

        if (vox.w != -1)
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
