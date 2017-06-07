#include "gpu_utils.cuh"

#include "gpu_errorcheck.cuh"

///////////////////////////////////////
template <class T>
void allocate_scorer(T*& s, size_t n)
{
    gpuErrchk( cudaMalloc( (void **) &s, sizeof(T)*n) );
    gpuErrchk( cudaMemset( (void *) s, 0, sizeof(T)*n) );
}
template void allocate_scorer<float4>(float4*&, size_t);
template void allocate_scorer<short2>(short2*&, size_t);
template void allocate_scorer<float>(float*&, size_t);

///////////////////////////////////////
template <class S, class T>
void retrieve_scorer(S* h, T* d, size_t n)
{
    gpuErrchk( cudaMemcpy(h, d, sizeof(T)*n, cudaMemcpyDeviceToHost) );
}
template void retrieve_scorer<float, float4>(float*, float4*, size_t);

///////////////////////////////////////
template <class T>
void array_to_device(T*& dest, const T* src, size_t n)
{
    gpuErrchk( cudaMalloc((void **) &dest, sizeof(T)*n) );
    gpuErrchk( cudaMemcpy(dest, src, sizeof(T)*n, cudaMemcpyHostToDevice) );
}
template void array_to_device<short>(short*&, const short*, size_t);
template void array_to_device<float2>(float2*&, const float2*, size_t);

///////////////////////////////////////

// template <class T>
// void symbol_to_device(T*& dest, const T* src, size_t n)
// {
//     gpuErrchk( cudaMalloc((void **) &dest, sizeof(T)*n) );
//     gpuErrchk( cudaMemcpyToSymbol(dest,  src, sizeof(T)*n, 0, cudaMemcpyHostToDevice) );
// }
// template void symbol_to_device<float4>(float4*& dest, const float4*, size_t);
// template void symbol_to_device<short2>(short2*& dest, const short2*, size_t);
