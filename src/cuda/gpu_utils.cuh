#ifndef __GPU_UTILS_CUH__
#define __GPU_UTILS_CUH__

#include "gpu_errorcheck.cuh"
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <vector>

template <class T>
void allocate_scorer(T*& s, size_t n);

template <class S, class T>
void retrieve_scorer(S* host, T* dev, size_t n);

template <class T, class S = T>
void array_to_device(T*& dest, const S* src, size_t n);

template <class T>
void symbol_to_device(T*& dest, const T* src, size_t n);

template <class T>
void sendVectorToTexture(size_t w, size_t h, size_t d,
                         std::vector<T> host_vec,
                         cudaArray* array,
                         texture<T, 3, cudaReadModeElementType>& tex);

template <class T>
__device__ __host__ void bubble_sort(T arr[], int n, int indx[]);

template <class T>
__device__ __host__ void swap(T* a, T* b);

__device__ void sum_kahan (float& old, const float& val, float& err);
__device__ void sum_kahan (float3& old, const float3& val, float3& err);
__device__ void sum_mul_kahan (float& old, const float& mult,
                               const float& val, float& err);
__device__ void sum_mul_kahan (float3& old, const float3& mult,
                               const float& val, float3& err);

#endif