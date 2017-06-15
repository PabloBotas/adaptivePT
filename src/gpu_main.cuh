#ifndef __GPU_MAIN_CUH__
#define __GPU_MAIN_CUH__

#include "patient_parameters.hpp"
#include "volume.hpp"

#include <iostream>
#include <string>
#include <vector>

#include <helper_cuda.h>
#include <helper_math.h>

void initialize_device (cudaEvent_t& start);

void stop_device (cudaEvent_t& start);

void gpu_raytrace (const Patient_Parameters_t& pat,
                   std::vector< Vector4_t<float> >& endpoints,
                   std::string output_file,
                   const std::vector< Vector4_t<float> >& orig_endpoints = std::vector< Vector4_t<float> >());

void gpu_raytrace_original (const Patient_Parameters_t &pat,
                            const Patient_Volume_t &ct,
                            std::vector< Vector4_t<float> >& endpoints,
                            std::vector< Vector4_t<float> >& init_pos,
                            std::vector< Vector4_t<float> >& init_pat_pos,
                            std::string output_file = std::string());

void gpu_raytrace_warped (const Patient_Parameters_t &pat,
                          const Patient_Volume_t &ct,
                          const std::vector< Vector4_t<float> >& orig_endpoints,
                          const std::vector< Vector4_t<float> >& init_pos,
                          std::vector< Vector4_t<float> >& endpoints,
                          std::string output_file = std::string());

void printDevProp (const int device, bool verbose);

#endif
