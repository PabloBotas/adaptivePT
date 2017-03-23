#ifndef __GPU_MATERIAL_CUH__
#define __GPU_MATERIAL_CUH__

#include <vector>
#include <string>

std::vector<float> readDensityCorrect(std::string f);
float HU2dens(int huValue);
//      convert HU to dens, in g/cm^3

#endif //__LIBMATERIAL_CU__
