#ifndef __GPU_MAIN_CUH__
#define __GPU_MAIN_CUH__

#include "patient_parameters.hpp"
#include "volume.hpp"

#include <iostream>
#include <string>
#include <vector>

#include <helper_cuda.h>
#include <helper_math.h>

void initialize_device(cudaEvent_t& start, cudaEvent_t& stop);

void stop_device(cudaEvent_t& start, cudaEvent_t& stop);

std::vector<float4> gpu_get_beam_endpoints(const Patient_Parameters_t &pat,
                                           const Patient_Volume_t &ct);

void runCalculation(const Patient_Parameters_t &pat,
                    const Patient_Volume_t &ct,
                    std::vector<float4>& restuls);

void printDevProp(const int device, bool verbose);

#endif
