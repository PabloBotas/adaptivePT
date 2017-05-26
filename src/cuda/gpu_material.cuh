#ifndef __GPU_MATERIAL_CUH__
#define __GPU_MATERIAL_CUH__

#include <vector>
#include <string>

std::vector<float> readDensityCorrect(const std::string f);
float HU2dens(const short val);
//      convert HU to dens, in g/cm^3
int HU2matId(const int val, const std::vector<int> hu_indexes);

#endif //__LIBMATERIAL_CU__
