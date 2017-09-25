#include "gpu_device_globals.cuh"

__device__ double4 xdata[NRAYS];  // x y z wepl (in cm)
__device__ double4 vxdata[NRAYS]; // vx vy vz (normalized) energy (eV)
__device__ short2 ixdata[NRAYS]; // particle metadata (beam and spot ID)

// CT variables
__device__ __constant__ double3 ctVoxSize;
__device__ __constant__ int3 ctVox;
__device__ __constant__ size_t ctTotalVoxN;

// Density-filled array
cudaArray* dens;
texture<float, 3, cudaReadModeElementType> dens_tex;

// Material ID-filled array
cudaArray* matid;
texture<float, 3, cudaReadModeElementType> matid_tex;

// Water resticted stopping power
__device__ __constant__ float stp_w_min_e;
__device__ __constant__ float stp_w_delta_e;
cudaArray *stp_w_array, *stp_w_b_coeff_array;
texture<float,1,cudaReadModeElementType> stp_w_tex, stp_w_b_coeff_tex;

// Stopping power ratio table
__device__ __constant__ float stp_ratio_min_e;
__device__ __constant__ float stp_ratio_delta_e;
cudaArray* stp_ratio_array;
texture<float, 2, cudaReadModeElementType> stp_ratio_tex;

// Bragg peaks LUTs
__device__ __constant__ float bp_energy_min;
__device__ __constant__ float bp_energy_delta;
__device__ __constant__ float bp_depth_delta;
cudaArray* bp_n_array;
texture<float, 2, cudaReadModeElementType> bp_n_tex;
cudaArray* bp_w_array;
texture<float, 2, cudaReadModeElementType> bp_w_tex;
cudaArray* bp_s_array;
texture<float, 2, cudaReadModeElementType> bp_s_tex;
cudaArray* bp_b_array;
texture<float, 2, cudaReadModeElementType> bp_b_tex;