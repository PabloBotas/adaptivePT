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
                   Array4<double>& endpoints,
                   std::string output_file,
                   const Array4<double>& orig_endpoints = Array4<double>());

void gpu_raytrace_original (const Patient_Parameters_t &pat,
                            const Volume_t &ct,
                            Array4<double>& endpoints,
                            Array4<double>& init_pos,
                            Array4<double>& init_pat_pos,
                            std::string output_file,
                            Array4<double>& influence);

void gpu_raytrace_warped (const Patient_Parameters_t &pat,
                          const Volume_t &ct,
                          const Array4<double>& orig_endpoints,
                          const Array4<double>& init_pos,
                          Array4<double>& endpoints,
                          std::string output_file,
                          Array4<double>& influence);

void gpu_calculate_influence (const short& pat,
                              const Array4<double>& endpoints,
                              Array4<double>& influence,
                              std::vector<float>& spot_weights,
                              std::vector<float>& inf_volume);

void printDevProp (const int device, bool verbose);

#endif
