#include "device_interaction.cuh"

//==========================================================
//      global variables
//==========================================================
// Rays
// device
__device__ float3 posbuffer[NRAYS]; // x y z    (in cm)
__device__ float3 dirbuffer[NRAYS]; // vx vy vz (normalized)
__device__ unsigned int mask_idx[NRAYS]; // vox numbers

// scoring array
float *scorer;

// CT variables
// device
cudaArray *dens;
texture<float,3,cudaReadModeElementType> dens_tex;
__device__ __constant__ float3 ctVoxSize;
__device__ __constant__ int3 ctVox;
__device__ __constant__ size_t ctTotalVoxN;
