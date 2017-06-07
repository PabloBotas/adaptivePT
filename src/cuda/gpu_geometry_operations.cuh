#ifndef __GPU_GEOMETRY_CUH__
#define __GPU_GEOMETRY_CUH__

#include <string>

#include "gpu_ray_class.cuh"

enum VoxelUpdater { UPDATEX, UPDATEY, UPDATEZ };
enum VoxelStepper { FORWARD = 1, BACKWARD = -1 };

__device__ __host__ float3 ext_to_int_coordinates(float3 a);
__device__ __host__ float4 ext_to_int_coordinates(float4 a);
__device__ __host__ float3 int_to_ext_coordinates(float3 a);
__device__ __host__ float4 int_to_ext_coordinates(float4 a);
__device__ float4 rotate(const float4& p, const float& gantry, const float& couch);
__device__ __host__ int getabs(int xvox, int yvox, int zvox, int ny, int nz);
__device__ float inters(const Ray& ray, const int4& vox, VoxelUpdater &indexvox, VoxelStepper &dvox);
__device__ void changeVoxel(int4 &vox, const VoxelUpdater indexvox, const VoxelStepper dvox);
__device__ int getVoxelIndex(int4 vox);
__device__ int getVoxelIndex(int3 vox);
__device__ int3 getVoxelCoords(unsigned int index);
__device__ float3 getVoxelCenter(int3 vox);
__device__ float3 getVoxelCenter(unsigned int index);

#endif  // GEOMETRY_CUH
