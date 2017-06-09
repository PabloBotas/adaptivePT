#ifndef __GPU_UTILS_CUH__
#define __GPU_UTILS_CUH__

#include "gpu_errorcheck.cuh"
#include <vector>

template <class T>
void allocate_scorer(T*& s, size_t n);

template <class S, class T>
void retrieve_scorer(S* host, T* dev, size_t n);

template <class T>
void array_to_device(T*& dest, const T* src, size_t n);

template <class T>
void symbol_to_device(T*& dest, const T* src, size_t n);

template <class T>
void sendVectorToTexture(size_t w, size_t h, size_t d,
                         std::vector<T> host_vec,
                         cudaArray* array,
                         texture<T, 3, cudaReadModeElementType>& tex);

#endif