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

std::vector< Vector4_t<float> > gpu_raytrace_plan(const Patient_Parameters_t &pat,
                                                  const Patient_Volume_t &ct);

void gpu_raytrace_plan(const Patient_Parameters_t &pat,
                       const Patient_Volume_t &ct,
                       std::vector< Vector4_t<float> >& restuls);

void printDevProp(const int device, bool verbose);

#endif
