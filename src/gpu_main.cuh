#ifndef __GPU_MAIN_CUH__
#define __GPU_MAIN_CUH__

#include "patient_parameters.hpp"
#include "command_line_parser.hpp"
#include "volume.hpp"

#include <iostream>
#include <string>
#include <vector>

#include <helper_cuda.h>
#include <helper_math.h>

void initialize_device (cudaEvent_t& start);

void stop_device (cudaEvent_t& start);

void gpu_raytrace (const Patient_Parameters_t& pat,
                   Array4<float>& endpoints,
                   std::string output_file,
                   const Array4<float>& orig_endpoints = Array4<float>());

void gpu_raytrace_original (const Patient_Parameters_t& pat,
                            const Parser& parser,
                            const Volume_t& ct,
                            Array4<float>& endpoints,
                            Array4<float>& initpos_xbuffer_dbg,
                            Array4<float>& initpos);

void gpu_raytrace_warped (const Patient_Parameters_t &pat,
                          const Parser& parser,
                          const Volume_t &ct,
                          const Array4<float>& orig_endpoints,
                          const Array4<float>& init_pos,
                          Array4<float>& endpoints,
                          std::vector<float>& new_energies);

void gpu_calculate_influence (const Patient_Parameters_t& pat,
                              const Parser& parser,
                              const Volume_metadata_t& ct_metadata,
                              Array4<float>& influence,
                              const std::vector<float>& new_energies = std::vector<float>());

void printDevProp (const int device, bool verbose);

#endif
