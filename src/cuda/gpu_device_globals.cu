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
texture<float, cudaTextureType3D, cudaReadModeElementType> dens_tex;

// Material ID-filled array
cudaArray* matid;
texture<float, cudaTextureType3D, cudaReadModeElementType> matid_tex;

// Boolean mask array
__device__ bool masking_vf;
cudaArray* vf_mask;
texture<int, cudaTextureType3D, cudaReadModeElementType> vf_mask_tex;

// Water resticted stopping power
__device__ __constant__ float stp_w_min_e;
__device__ __constant__ float stp_w_delta_e;
cudaArray *stp_w_array, *stp_w_b_coeff_array;
texture<float, cudaTextureType1D, cudaReadModeElementType> stp_w_tex, stp_w_b_coeff_tex;

// Stopping power ratio table
__device__ __constant__ float stp_ratio_min_e;
__device__ __constant__ float stp_ratio_delta_e;
cudaArray* stp_ratio_array;
texture<float, cudaTextureType2D, cudaReadModeElementType> stp_ratio_tex;

// Bragg peaks LUTs
__device__ __constant__ float bp_energy_min;
__device__ __constant__ float bp_energy_delta;
__device__ __constant__ float bp_depth_delta;
cudaArray* bp_n_array;
texture<float, cudaTextureType2D, cudaReadModeElementType> bp_n_tex;
cudaArray* bp_w_array;
texture<float, cudaTextureType2D, cudaReadModeElementType> bp_w_tex;
cudaArray* bp_s_array;
texture<float, cudaTextureType2D, cudaReadModeElementType> bp_s_tex;
cudaArray* bp_b_array;
texture<float, cudaTextureType2D, cudaReadModeElementType> bp_b_tex;
cudaArray* bp_range_array;
texture<float, cudaTextureType1D, cudaReadModeElementType> bp_range_tex;
