#include "gpu_utils.cuh"

#include "gpu_errorcheck.cuh"
#include "vector4.hpp"

///////////////////////////////////////
template <class T>
void allocate_scorer(T*& s, size_t n)
{
    gpuErrchk( cudaMalloc( (void **) &s, sizeof(T)*n) );
    gpuErrchk( cudaMemset( (void *) s, 0, sizeof(T)*n) );
}
template void allocate_scorer<float4>(float4*&, size_t);
template void allocate_scorer<float3>(float3*&, size_t);
template void allocate_scorer<short2>(short2*&, size_t);
template void allocate_scorer<float>(float*&, size_t);
template void allocate_scorer<unsigned long long int>(unsigned long long int*&, size_t);

///////////////////////////////////////
template <class S, class T>
void retrieve_scorer(S* h, T* d, size_t n)
{
    gpuErrchk( cudaMemcpy(h, d, sizeof(T)*n, cudaMemcpyDeviceToHost) );
}
template void retrieve_scorer<short, short2>(short*, short2*, size_t);
template void retrieve_scorer<float, float>(float*, float*, size_t);
template void retrieve_scorer<float, float4>(float*, float4*, size_t);
template void retrieve_scorer<float, float3>(float*, float3*, size_t);
template void retrieve_scorer<unsigned long long int, unsigned long long int>(unsigned long long int*, unsigned long long int*, size_t);

///////////////////////////////////////
template <class T, class S>
void array_to_device(T*& dest, const S* src, size_t n)
{
    gpuErrchk( cudaMalloc((void **) &dest, sizeof(T)*n) );
    gpuErrchk( cudaMemcpy(dest, src, sizeof(T)*n, cudaMemcpyHostToDevice) );
}
template void array_to_device<short>(short*&, const short*, size_t);
template void array_to_device<float>(float*&, const float*, size_t);
template void array_to_device<float2>(float2*&, const float2*, size_t);
template void array_to_device<float3>(float3*&, const float3*, size_t);
template void array_to_device<float4>(float4*&, const float4*, size_t);
template void array_to_device<float4, Vector4_t<float>>(float4*&, const Vector4_t<float>*, size_t);
template void array_to_device<float3, Vector3_t<float>>(float3*&, const Vector3_t<float>*, size_t);

///////////////////////////////////////

// template <class T>
// void symbol_to_device(T*& dest, const T* src, size_t n)
// {
//     gpuErrchk( cudaMalloc((void **) &dest, sizeof(T)*n) );
//     gpuErrchk( cudaMemcpyToSymbol(dest,  src, sizeof(T)*n, 0, cudaMemcpyHostToDevice) );
// }
// template void symbol_to_device<float4>(float4*& dest, const float4*, size_t);
// template void symbol_to_device<short2>(short2*& dest, const short2*, size_t);

///////////////////////////////////////
template <class T>
void sendVectorToTexture(size_t w, size_t h, size_t d,
                         std::vector<T> host_vec,
                         cudaArray* array,
                         texture<T, 3, cudaReadModeElementType>& tex)
{
    //  create a 3d array on device
    cudaExtent extent = make_cudaExtent(w, h, d);
    gpuErrchk( cudaMalloc3DArray(&array, &tex.channelDesc, extent) );

    // copy data to GPU
    cudaMemcpy3DParms pars = cudaMemcpy3DParms();
    pars.srcPtr   = make_cudaPitchedPtr((void *)host_vec.data(),
                                        extent.width*sizeof(T),
                                        extent.width, extent.height);
    pars.dstArray = array;
    pars.extent   = extent;
    pars.kind     = cudaMemcpyHostToDevice;
    gpuErrchk( cudaMemcpy3D(&pars) );
    // Bind device array to texture
    tex.normalized = false;
    tex.filterMode = cudaFilterModePoint;
    gpuErrchk( cudaBindTextureToArray(tex, array, tex.channelDesc) );
}
template void sendVectorToTexture<float>(size_t w, size_t h, size_t d,
                                         std::vector<float> host_vec,
                                         cudaArray* array,
                                         texture<float, 3, cudaReadModeElementType>& tex);
template void sendVectorToTexture<int>(size_t w, size_t h, size_t d,
                                       std::vector<int> host_vec,
                                       cudaArray* array,
                                       texture<int, 3, cudaReadModeElementType>& tex);
///////////////////////////////////////
template <class T>
__device__ __host__ void bubble_sort(T arr[], int n, int indx[])
{
    for (int i = 0; i < n-1; i++)      
        for (int j = 0; j < n-i-1; j++) 
            if (arr[j] > arr[j+1]) {
                swap<T>(&arr[j], &arr[j+1]);
                if (indx != NULL)
                    swap<int>(&indx[j], &indx[j+1]);
            }
}
template void bubble_sort(float arr[], int n, int indx[]);

template <class T>
__device__ __host__ void swap(T* a, T* b)
{
    T t = *a;
    *a = *b;
    *b = t;
}
template void swap(int*, int*);
template void swap(float*, float*);

///////////////////////////////////////
__device__ void sum_kahan (float& old, const float& val, float& err) {
    float change = val - err;
    float new_val = old + val;
    err = (new_val - old) - change;
    old = new_val;
}

__device__ void sum_kahan (float3& old, const float3& val, float3& err) {
    sum_kahan(old.x, val.x, err.x);
    sum_kahan(old.y, val.y, err.y);
    sum_kahan(old.z, val.z, err.z);
}

__device__ void sum_mul_kahan (float& old, const float& mult,
                               const float& val, float& err) {
    float change = __fmaf_rn(mult, val, -err);
    float new_val = __fmaf_rn(mult, val, old);
    err = (new_val - old) - change;
    old = new_val;
}

__device__ void sum_mul_kahan (float3& old, const float3& mult,
                               const float& val, float3& err) {
    sum_mul_kahan(old.x, mult.x, val, err.x);
    sum_mul_kahan(old.y, mult.y, val, err.y);
    sum_mul_kahan(old.z, mult.z, val, err.z);
}
