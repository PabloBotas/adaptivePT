#ifndef __DEVICE_INTERACTION_CUH__
#define __DEVICE_INTERACTION_CUH__

// includes, project
#include <cuda_runtime.h>

// helper
#include <helper_cuda.h>
#include <helper_math.h>

//==========================================================
//      GPU configurations
//==========================================================
//      the leading dimension of the 2d thread grid: total registers per block
#define NBLOCKX 32768
//      number of threads per block
#define NTHREAD_PER_BLOCK_SOURCE 256 // set source direction
#define NTHREAD_PER_BLOCK_RAYS   64  // rays

//==========================================================
//      mathematical constants
//==========================================================
#define PI 3.1415926535897932384626433
#define ZERO 1.0e-20
#define SZERO 1.0e-6
#define INF 1.0e20
#define MM2CM 0.1
#define CM2MM 10
#define MeV2eV 1000000

//==========================================================
//      global variables
//==========================================================
#define NRAYS 32768 // max number of rays calculated simultaneously

extern __device__ float4 xdata[NRAYS];  // x y z wepl (in cm)
extern __device__ float4 vxdata[NRAYS]; // vx vy vz (normalized) energy (eV)
extern __device__ short2 ixdata[NRAYS]; // particle metadata (beam and spot ID)

// CT variables
extern __device__ __constant__ float3 ctVoxSize;
extern __device__ __constant__ int3 ctVox;
extern __device__ __constant__ size_t ctTotalVoxN;

// Density-filled array
extern cudaArray* dens;
extern texture<float, 3, cudaReadModeElementType> dens_tex;

// Material ID-filled array
extern cudaArray* matid;
extern texture<float, 3, cudaReadModeElementType> matid_tex;

// Stopping power table
extern __device__ __constant__ float stp_ratio_min_e;
extern __device__ __constant__ float stp_ratio_delta_e;
extern cudaArray* stp_ratio_array;
extern texture<float, 2, cudaReadModeElementType> stp_ratio_tex;

#endif
