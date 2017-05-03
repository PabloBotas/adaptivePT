#ifndef __PATIENT_VOLUME_HPP__
#define __PATIENT_VOLUME_HPP__

#include "special_types.hpp"

#include <string>
#include <vector>

class Patient_Volume_t
{
public:
    enum Source_type
    {
        CTVOLUME,
        MHA
    };

    Patient_Volume_t(std::string file,
                    Patient_Volume_t::Source_type source);
    Patient_Volume_t(std::string file,
                     Patient_Volume_t::Source_type source,
                     unsigned int nx, unsigned int ny, unsigned int nz,
                     float dx, float dy, float dz);
    Patient_Volume_t(const float* src,
    		         const CT_Dims_t dims);
    Patient_Volume_t(const CT_Dims_t dims);

    void setVoxels(unsigned int x, unsigned int y, unsigned int z);
    void setSpacing(float x, float y, float z);
    void output(std::string outfile);

    Patient_Volume_t::Source_type source_type;
    std::string file;
    Vector_t<unsigned int> n;
    Vector_t<float> d;
    Vector_t<float> origin;
    Vector_t<float> imgCenter;
    unsigned int nElements;

    std::vector<float> hu;
private:
    void read_volume();
};

#endif
