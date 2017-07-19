#include "gpu_geometry_tools.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_ray_class.cuh"

#include "vector_types.h"
#include <assert.h>

// CT Navigation ----------------------------
__device__ double to_boundary(const double3& pos,
                              const double3& dir,
                              const int4& vox,
                              VoxelUpdater& voxUpdater,
                              VoxelStepper& voxStepper)
{
    //    Checking out all the voxel walls for the smallest distance...
    // Z
    double invcos = (dir.z != 0.0f) ? 1.0f/dir.z : INF;
    int ifNext = (invcos > 0.0) ? 1 : 0;
    double step = ((vox.z+ifNext)*ctVoxSize.z - pos.z) * invcos;
    voxUpdater = UPDATEZ;
    voxStepper = ifNext ? FORWARD : BACKWARD;
    // Y
    invcos = (dir.y != 0.0f) ? 1.0f/dir.y : INF;
    ifNext = (invcos > 0.0f) ? 1 : 0;
    double tempstep = ((vox.y+ifNext) * ctVoxSize.y - pos.y) * invcos;
    if (tempstep < step)
    {
        step = tempstep;
        voxUpdater = UPDATEY;
        voxStepper = ifNext ? FORWARD : BACKWARD;
    }
    // X
    invcos = (dir.x != 0.0f) ? 1.0f/dir.x : INF;
    ifNext = (invcos > 0.0f) ? 1 : 0;
    tempstep = ((vox.x+ifNext) * ctVoxSize.x - pos.x) * invcos;
    if (tempstep < step)
    {
        step = tempstep;
        voxUpdater = UPDATEX;
        voxStepper = ifNext ? FORWARD : BACKWARD;
    }
    return fabs(step);
    // int3 v = make_int3(vox);
    // double3 invcos = inverse(dir);
    // int3 ifNext = bigger(invcos, 0.0);
    // double3 step = ((v+ifNext)*ctVoxSize - pos) * invcos;
    // int i = find_min(step);
    // voxUpdater = static_cast<VoxelUpdater>(i);
    // voxStepper = static_cast<VoxelStepper>(at(ifNext, i));
    // return(fabs(fminf(step)));
}

__device__ double to_boundary(const double3& pos,
                              const double3& dir,
                              const int4& vox,
                              VoxelUpdater& voxUpdater,
                              VoxelStepper& voxStepper,
                              const double3& endpoint)
{
    double boundary = to_boundary(pos, dir, vox, voxUpdater, voxStepper);
    double const min_dist = 0.00003;

    double3 r = endpoint-pos;
    double dist = length(r);
    double cos_to_point = dot(r, dir)/(dist*length(dir));

    if(dist > min_dist &&
       !(cos_to_point >  0.9999 && cos_to_point < 1.0001) &&
       !(cos_to_point > -0.0001 && cos_to_point < 0.0001))
    {
        int i =blockIdx.x*blockDim.x + threadIdx.x;
        printf("WARNING! %d - %f %f - %f %f %f - %f %f %f - %f %f %f - %f %f %f - %f %f %f\n", 
               i, cos_to_point, dist,
               xdata[i].x, xdata[i].y, xdata[i].z,
               vxdata[i].x, vxdata[i].y, vxdata[i].z,
               pos.x, pos.y, pos.z,
               dir.x, dir.y, dir.z,
               endpoint.x, endpoint.y, endpoint.z);
        assert((cos_to_point >  0.9999 && cos_to_point < 1.0001) ||
               (cos_to_point > -0.0001 && cos_to_point < 0.0001));
    }

    if(dist <= min_dist)
    {
        boundary = 0;
        voxUpdater = NONE;
    }
    else if(dist < boundary)
    {
        boundary = dist;
        voxUpdater = NONE;
    }

    // printf("Dist to endpoint - step - cos: %f - %f\n", dist, boundary);

    return boundary;
}

__device__ int ahead_or_behind(const double3& dir,
                               const double3& point,
                               const double3& pos)
{
    double3 vector = point - pos;
    return dot(vector, dir) >= 0 ? 1 : -1;
}

__device__ void changeVoxel(int4& vox,
                            const VoxelUpdater updater,
                            const VoxelStepper stepper)
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
    else if (updater == UPDATEX)
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

__device__ int4 get_voxel (double3 pos)
{
    int4 vox;
    vox.x = floor(pos.x/ctVoxSize.x);
    vox.y = floor(pos.y/ctVoxSize.y);
    vox.z = floor(pos.z/ctVoxSize.z);
    // Check if in CT grid
    if (vox.x < 0 || vox.x >= ctVox.x ||
        vox.y < 0 || vox.y >= ctVox.y ||
        vox.z < 0 || vox.z >= ctVox.z)
        vox.w = -1;
    else
        vox.w = vox.z + vox.y*ctVox.z + vox.x*ctVox.z*ctVox.y;

    return vox;
}

__device__ int3 getVoxelCoords(unsigned int index)
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


__device__ double3 getVoxelCenter(int3 vox)
{
    double3 pos;
    pos.x = vox.x*ctVoxSize.x + 0.5*ctVoxSize.x;
    pos.y = vox.y*ctVoxSize.y + 0.5*ctVoxSize.y;
    pos.z = vox.z*ctVoxSize.z + 0.5*ctVoxSize.z;

    return pos;
}

__device__ double3 getVoxelCenter(unsigned int index)
{
    return getVoxelCenter(getVoxelCoords(index));
}

