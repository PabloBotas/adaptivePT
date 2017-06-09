#ifndef __GPU_MAIN_CUH__
#define __GPU_MAIN_CUH__

#include "patient_parameters.hpp"
#include "volume.hpp"

#include <iostream>
#include <string>
#include <vector>

#include <helper_cuda.h>
#include <helper_math.h>

void initialize_device(cudaEvent_t& start);

void stop_device(cudaEvent_t& start);

void gpu_raytrace_plan(const Patient_Parameters_t &pat,
                       const Patient_Volume_t &ct,
                       std::vector< Vector4_t<float> >& endpoints,
                       std::vector< Vector4_t<float> >& directions,
                       std::vector< Vector2_t<short> >& metadata);

void
gpu_backtrace_endpoints(const Patient_Parameters_t& pat,
                        const Patient_Volume_t& ct,
                        const std::vector< Vector4_t<float> >& endpoints,
                        const std::vector< Vector4_t<float> >& directions,
                        const std::vector< Vector2_t<short> >& metadata,
                        std::vector< Vector4_t<float> >& output);

void printDevProp(const int device, bool verbose);

#endif
