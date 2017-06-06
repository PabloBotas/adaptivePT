#ifndef __GPU_UTILS_CUH__
#define __GPU_UTILS_CUH__

#include "gpu_errorcheck.cuh"

template <class T>
void allocate_scorer(T*& s, size_t n);

template <class S, class T>
void retrieve_scorer(S* host, T* dev, size_t n);

#endif