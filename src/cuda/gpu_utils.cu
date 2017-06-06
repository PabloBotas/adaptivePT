#include "gpu_utils.cuh"

#include "gpu_errorcheck.cuh"

template <class T>
void allocate_scorer(T*& s, size_t n)
{
    gpuErrchk( cudaMalloc( (void **) &s, sizeof(T)*n) );
    gpuErrchk( cudaMemset( (void *) s, 0, sizeof(T)*n) );
}

template void allocate_scorer<float4>(float4*&, size_t);
template void allocate_scorer<short2>(short2*&, size_t);
template void allocate_scorer<float>(float*&, size_t);

template <class S, class T>
void retrieve_scorer(S* h, T* d, size_t n)
{
    gpuErrchk( cudaMemcpy(h, d, sizeof(T)*n, cudaMemcpyDeviceToHost) );
}

template void retrieve_scorer<float, float4>(float*, float4*, size_t);

// template void retrieve_scorer<short, short2>(short*& host, short2* dev, size_t n);
// template void retrieve_scorer<float, float>(float* host, float* dev, size_t n);
