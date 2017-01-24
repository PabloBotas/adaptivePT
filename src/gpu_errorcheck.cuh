#ifndef __GPU_ERRORCHECK_CUH__
#define __GPU_ERRORCHECK_CUH__

#include <string>

#define gpuErrchk(ans) {gpuAssert((ans), __FILE__, __LINE__); }

void gpuAssert(cudaError_t code, char const *file, int line, bool abort=true);

void ioError(std::string info);

#endif
