#include "gpu_geometry_tools.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_ray_class.cuh"

#include "vector_types.h"
#include <assert.h>

// CT Navigation ----------------------------
__device__ float to_boundary(const float3& pos,
                             const float3& dir,
                             const int4& vox,
                             VoxelUpdater& voxUpdater,
                             VoxelStepper& voxStepper)
{
    //    Checking out all the voxel walls for the smallest distance...
    // Z
    float invcos = (dir.z != 0.0f) ? 1.0f/dir.z : INF;
    int ifNext = (invcos > 0.0) ? 1 : 0;
    float step = ((vox.z+ifNext)*ctVoxSize.z - pos.z) * invcos;
    voxUpdater = UPDATEZ;
    voxStepper = ifNext ? FORWARD : BACKWARD;
    // Y
    invcos = (dir.y != 0.0f) ? 1.0f/dir.y : INF;
    ifNext = (invcos > 0.0f) ? 1 : 0;
    float tempstep = ((vox.y+ifNext) * ctVoxSize.y - pos.y) * invcos;
    if (tempstep < step) {
        step = tempstep;
        voxUpdater = UPDATEY;
        voxStepper = ifNext ? FORWARD : BACKWARD;
    }
    // X
    invcos = (dir.x != 0.0f) ? 1.0f/dir.x : INF;
    ifNext = (invcos > 0.0f) ? 1 : 0;
    tempstep = ((vox.x+ifNext) * ctVoxSize.x - pos.x) * invcos;
    if (tempstep < step) {
        step = tempstep;
        voxUpdater = UPDATEX;
        voxStepper = ifNext ? FORWARD : BACKWARD;
    }
    return fabs(step);
    // int3 v = make_int3(vox);
    // float3 invcos = inverse(dir);
    // int3 ifNext = bigger(invcos, 0.0);
    // float3 step = ((v+ifNext)*ctVoxSize - pos) * invcos;
    // int i = find_min(step);
    // voxUpdater = static_cast<VoxelUpdater>(i);
    // voxStepper = static_cast<VoxelStepper>(at(ifNext, i));
    // return(fabs(fminf(step)));
}

__device__ float to_boundary(const float3& pos,
                             const float3& dir,
                             const int4& vox,
                             VoxelUpdater& voxUpdater,
                             VoxelStepper& voxStepper,
                             const float3& endpoint)
{
    float boundary = to_boundary(pos, dir, vox, voxUpdater, voxStepper);

    float3 r = endpoint-pos;
    float dist = length(r);
    float cos_to_point = dot(r, dir)/(dist*length(dir));

    float const min_dist = 0.005;
    if (dist > min_dist &&
       !(cos_to_point >  0.99 && cos_to_point < 1.01) &&
       !(cos_to_point > -0.01 && cos_to_point < 0.01)) {
        int i = blockIdx.x*blockDim.x + threadIdx.x;
        printf("WARNING! %d - %f %f - %f %f %f - %f %f %f - %f %f %f - %f %f %f - %f %f %f\n", 
               i, cos_to_point, dist,
               xdata[i].x, xdata[i].y, xdata[i].z,
               vxdata[i].x, vxdata[i].y, vxdata[i].z,
               pos.x, pos.y, pos.z,
               dir.x, dir.y, dir.z,
               endpoint.x, endpoint.y, endpoint.z);
        assert((cos_to_point >  0.9995 && cos_to_point < 1.0005) ||
               (cos_to_point > -0.0005 && cos_to_point < 0.0005));
    }

    if (dist <= min_dist) {
        boundary = 0;
        voxUpdater = NONE;
    }
    else if(dist < boundary) {
        boundary = dist;
        voxUpdater = NONE;
    }

    // printf("Dist to endpoint - step - cos: %f - %f\n", dist, boundary);

    return boundary;
}

__device__ int ahead_or_behind(const float3& dir,
                               const float3& point,
                               const float3& pos)
{
    float3 vector = point - pos;
    return dot(vector, dir) >= 0 ? 1 : -1;
}

__device__ void changeVoxel(int4& vox,
                            const VoxelUpdater updater,
                            const VoxelStepper stepper)
//    Changes voxel according to the information passed by inters()
{
    if (updater == UPDATEZ) {
        vox.z += stepper;
        vox.w += stepper;
    }
    else if (updater == UPDATEY) {
        vox.y += stepper;
        vox.w += stepper*ctVox.z;
    }
    else if (updater == UPDATEX) {
        vox.x += stepper;
        vox.w += stepper*ctVox.z*ctVox.y;
    }

    if (vox.x < 0 || vox.x > ctVox.x ||
        vox.y < 0 || vox.y > ctVox.y ||
        vox.z < 0 || vox.z > ctVox.z)
        vox.w = -1;
}

__host__ __device__ int getabs(int xvox, int yvox, int zvox, int ny, int nz)
//    Gets the absolute voxel # from the coordinate voxel #s
{
    return zvox + yvox*nz + xvox*nz*ny;
}

__device__ int4 get_voxel (float3 pos)
{
    int4 vox;
    vox.x = int(floor(pos.x/ctVoxSize.x)+0.5);
    vox.y = int(floor(pos.y/ctVoxSize.y)+0.5);
    vox.z = int(floor(pos.z/ctVoxSize.z)+0.5);
    // Check if in CT grid
    if (vox.x < 0 || vox.x >= ctVox.x ||
        vox.y < 0 || vox.y >= ctVox.y ||
        vox.z < 0 || vox.z >= ctVox.z)
        vox.w = -1;
    else
        vox.w = vox.z + vox.y*ctVox.z + vox.x*ctVox.z*ctVox.y;

    return vox;
}


__device__ int get_voxel_abs (float3 pos)
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

    return vox.w;
}


__device__ int3 getVoxelCoords(unsigned int index)
{
    int3 vox = make_int3(-1,-1,-1);
    if (index < ctTotalVoxN) {
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

