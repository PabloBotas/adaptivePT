#include "gpu_device_globals.cuh"

__device__ float4 xdata[NRAYS];  // x y z wepl (in cm)
__device__ float4 vxdata[NRAYS]; // vx vy vz (normalized) energy (eV)
__device__ short2 ixdata[NRAYS]; // particle metadata (beam and spot ID)

// CT variables
__device__ __constant__ float3 ctVoxSize;
__device__ __constant__ int3 ctVox;
__device__ __constant__ size_t ctTotalVoxN;

// Density-filled array
cudaArray* dens;
texture<float, 3, cudaReadModeElementType> dens_tex;

// Material ID-filled array
cudaArray* matid;
texture<float, 3, cudaReadModeElementType> matid_tex;

// Stopping power table
__device__ __constant__ float stp_ratio_min_e;
__device__ __constant__ float stp_ratio_delta_e;
cudaArray* stp_ratio_array;
texture<float, 2, cudaReadModeElementType> stp_ratio_tex;