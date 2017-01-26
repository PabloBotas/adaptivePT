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

    Patient_Volume_t();
    Patient_Volume_t(std::string file,
                    Patient_Volume_t::Source_type source);
    Patient_Volume_t(std::string file,
                     Patient_Volume_t::Source_type source,
                     unsigned int nx, unsigned int ny, unsigned int nz,
                     float dx, float dy, float dz);

    void setVoxels(unsigned int x, unsigned int y, unsigned int z);
    void setSpacing(float x, float y, float z);

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