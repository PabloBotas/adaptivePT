#include "ray_kernel.cuh"
#include "device_interaction.cuh"
#include "ray_class.cuh"
#include "geometry_operations.cuh"

__global__ void calculateRays_kernel(int num, float *scorer)
{
    // const int id = blockIdx.x*blockDim.x + threadIdx.x;
    // if(id < num)
    // {
    //     Ray ray(posbuffer[id], dirbuffer[id]);
    //     int4 vox   = make_int4(floorf(ray.pos/ctVoxSize));
    //     vox.w = getVoxelIndex(vox);

    //     float wepl = 0.f;
    //     VoxelUpdater voxUpdater;
    //     VoxelStepper voxStepper;
    //     while (true)
    //     {
    //         float density = tex3D(dens_tex, vox.z, vox.y, vox.x);
    //         float step    = inters(ray, vox, voxUpdater, voxStepper);

    //         wepl += density*step;
    //         if(mask_idx[id] == vox.w)
    //         {
    //             scorer[vox.w] = wepl;
    //             break;
    //         }

    //         ray.step(step);
    //         changeVoxel(vox, voxUpdater, voxStepper);
    //     }
    // }
}
