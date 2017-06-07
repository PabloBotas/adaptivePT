#ifndef __GPU_GEOMETRY_CUH__
#define __GPU_GEOMETRY_CUH__

#include <string>

#include "gpu_ray_class.cuh"

enum VoxelUpdater { UPDATEX, UPDATEY, UPDATEZ };
enum VoxelStepper { FORWARD = 1, BACKWARD = -1 };

// CT navigation
__device__ __host__ int getabs(int xvox, int yvox, int zvox, int ny, int nz);
__device__ float inters(const Ray& ray, const int4& vox, VoxelUpdater &indexvox, VoxelStepper &dvox);
__device__ void changeVoxel(int4 &vox, const VoxelUpdater indexvox, const VoxelStepper dvox);
__device__ int getVoxelIndex(int4 vox);
__device__ int getVoxelIndex(int3 vox);
__device__ int3 getVoxelCoords(unsigned int index);
__device__ float3 getVoxelCenter(int3 vox);
__device__ float3 getVoxelCenter(unsigned int index);

#endif  // GEOMETRY_CUH
