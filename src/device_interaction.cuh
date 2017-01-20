#ifndef __DEVICE_INTERACTION_CUH__
#define __DEVICE_INTERACTION_CUH__

// includes, system
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <sys/stat.h>

// includes, project
#include <cuda_runtime.h>

// helper
#include <helper_cuda.h>
#include <helper_math.h>

//include stl
#include <vector>
#include <string>

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

//==========================================================
//      global variables
//==========================================================

#define NRAYS 131072 // number of rays calculated simultaneously

struct inputParFile_t {
    std::string densCorrFile;
    std::string mask;
    // give default values
    inputParFile_t()
    {
        mask = "CTStatic/mask.dat";
    }
};

#ifdef SEPARATE_COMPILATION

// Rays
// host
extern float3 posInit_h;
extern float3 dirInit_h;
// device
extern __device__ float3 posbuffer[NRAYS]; // x y z    (in cm)
extern __device__ float3 dirbuffer[NRAYS]; // vx vy vz (normalized)
extern __device__ unsigned int mask_idx[NRAYS]; // vox numbers

// scoring array
extern float *scorer;

// CT variables
// host
extern std::vector<short> ctnum_h;
extern std::vector<float> dens_h;
extern float3 ctVoxSize_h;
extern int3 ctVox_h;
extern size_t ctTotalVoxN_h = 0;
// device
extern cudaArray *dens;
extern texture<float,3,cudaReadModeElementType> dens_tex;
extern __device__ __constant__ float3 ctVoxSize;
extern __device__ __constant__ int3 ctVox;
extern __device__ __constant__ size_t ctTotalVoxN;

// Treatment variables
extern int nbeams_h;
extern std::vector<string> beamfolders_h;
extern std::vector<string> beamnames_h;
extern std::vector<float2> beamAngles_h;
extern std::vector<float>  isoToBeam_h;

#endif

#endif // GLOBAL_CUH
