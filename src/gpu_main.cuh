#ifndef __GPU_MAIN_CUH__
#define __GPU_MAIN_CUH__

#include "patient_parameters.hpp"
#include "volume.hpp"

#include <iostream>
#include <string>
#include <vector>


void gpu_main(const Patient_Parameters_t& pat, const Patient_Volume_t &ct);
void runCalculation(const Patient_Parameters_t &pat, const Patient_Volume_t &ct);
void printDevProp(const int device);

#endif