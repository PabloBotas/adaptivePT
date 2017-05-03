#include "gpu_device_interaction.cuh"

__device__ float4 xdata[NRAYS];  // x y z wepl (in cm)
__device__ float4 vxdata[NRAYS]; // vx vy vz (normalized) energy (eV)
__device__ short2 ixdata[NRAYS]; // particle metadata (beam and spot ID)

// CT variables
cudaArray *dens;
texture<float, 3, cudaReadModeElementType> dens_tex;
__device__ __constant__ float3 ctVoxSize;
__device__ __constant__ int3   ctVox;
__device__ __constant__ size_t ctTotalVoxN;
