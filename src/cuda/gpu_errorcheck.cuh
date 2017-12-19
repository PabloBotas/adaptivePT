#ifndef __GPU_ERRORCHECK_CUH__
#define __GPU_ERRORCHECK_CUH__

#include <string>
#include <cuda.h>
#include <cuda_runtime.h>

#define gpuErrchk(ans) {gpuAssert((ans), __FILE__, __LINE__); }

void gpuAssert(cudaError_t code, char const *file, int line, bool abort=true);
void check_kernel_execution(char const *file, int line);

#endif
