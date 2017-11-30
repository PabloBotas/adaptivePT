#include "gpu_errorcheck.cuh"
#include "cuda_helper/helper_cuda.h"

#include <iostream>
#include <string>

#define gpuErrchk(ans) {gpuAssert((ans), __FILE__, __LINE__); }

void gpuAssert(cudaError_t code, char const *file, int line, bool abort)
{
    if (code != cudaSuccess) {
        std::cerr << "Error code " << code << ": " << cudaGetErrorString(code) << std::endl;
        std::cerr << "File: " << file << std::endl;
        std::cerr << "Line: " << line << std::endl;
        if (abort)
            exit(code);
    }
}

void check_kernel_execution(char const *file, int line)
{
    cudaError_t code = cudaThreadSynchronize();
    gpuAssert(code, file, line);
    code = cudaPeekAtLastError();
    gpuAssert(code, file, line);
}
