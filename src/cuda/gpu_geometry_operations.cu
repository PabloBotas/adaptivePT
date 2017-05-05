#include "gpu_geometry_operations.cuh"
#include "gpu_device_interaction.cuh"
#include "gpu_ray_class.cuh"

#include "vector_types.h"

__device__ float inters(const Ray& ray,
                        const int4& vox,
                        VoxelUpdater& voxUpdater,
                        VoxelStepper& voxStepper)
{
    //    Checking out all the voxel walls for the smallest distance...
    // Z
    float invcos = (ray.dir.z != 0.0f) ? 1.0f/ray.dir.z : INF;
    int ifNext = (invcos > 0.0) ? 1 : 0;
    float step = ((vox.z+ifNext)*ctVoxSize.z - ray.pos.z) * invcos;
    voxUpdater = UPDATEZ;
    voxStepper = ifNext ? FORWARD : BACKWARD;
    // Y
    invcos = (ray.dir.y != 0.0f) ? 1.0f/ray.dir.y : INF;
    ifNext = (invcos > 0.0) ? 1 : 0;
    float tempstep = ((vox.y+ifNext) * ctVoxSize.y - ray.pos.y) * invcos;
    if (tempstep < step)
    {
        step = tempstep;
        voxUpdater = UPDATEY;
        voxStepper = ifNext ? FORWARD : BACKWARD;

    }
    // X
    invcos = (ray.dir.x != 0.0f) ? 1.0f/ray.dir.x : INF;
    ifNext = (invcos > 0.0) ? 1 : 0;
    tempstep = ((vox.x+ifNext) * ctVoxSize.x - ray.pos.x) * invcos;
    if (tempstep < step)
    {
        step = tempstep;
        voxUpdater = UPDATEX;
        voxStepper = ifNext ? FORWARD : BACKWARD;
    }

    return fabs(step);
}

__device__ void changeVoxel(int4& vox, const VoxelUpdater updater, const VoxelStepper stepper)
//    Changes voxel according to the information passed by inters()
{
    if (updater == UPDATEZ)
    {
        vox.z += stepper;
        vox.w += stepper;
    }
    else if (updater == UPDATEY)
    {
        vox.y += stepper;
        vox.w += stepper*ctVox.z;
    }
    else
    {
        vox.x += stepper;
        vox.w += stepper*ctVox.z*ctVox.y;
    }

    if(vox.x < 0 || vox.x > ctVox.x ||
       vox.y < 0 || vox.y > ctVox.y ||
       vox.z < 0 || vox.z > ctVox.z)
        vox.w = -1;
}

__host__ __device__ int getabs(int xvox, int yvox, int zvox, int ny, int nz)
//    Gets the absolute voxel # from the coordinate voxel #s
{
    return zvox + yvox*nz + xvox*nz*ny;
}

__device__ int getVoxelIndex(int3 vox)
//      return the absolute vox index according to the coordinate
{
    // Check if in CT grid
    if (vox.x < 0 || vox.x >= ctVox.x ||
        vox.y < 0 || vox.y >= ctVox.y ||
        vox.z < 0 || vox.z >= ctVox.z)
        return -1;

    return vox.z + vox.y*ctVox.z + vox.x*ctVox.z*ctVox.y;
}

__device__ int getVoxelIndex(int4 vox)
{
    return getVoxelIndex(make_int3(vox));
}

__device__ int3 getVoxelCoords(unsigned int index)
//      return the absolute vox index according to the coordinate
{
    int3 vox = make_int3(-1,-1,-1);
    if(index < ctTotalVoxN)
    {
        vox.x = index/(ctVox.y*ctVox.z) % ctVox.x;
        vox.y = index/ctVox.z % ctVox.y;
        vox.z = index % ctVox.z;
    }

    return vox;
}


__device__ float3 getVoxelCenter(int3 vox)
{
    float3 pos;
    pos.x = vox.x*ctVoxSize.x + 0.5*ctVoxSize.x;
    pos.y = vox.y*ctVoxSize.y + 0.5*ctVoxSize.y;
    pos.z = vox.z*ctVoxSize.z + 0.5*ctVoxSize.z;

    return pos;
}

__device__ float3 getVoxelCenter(unsigned int index)
{
    return getVoxelCenter(getVoxelCoords(index));
}

