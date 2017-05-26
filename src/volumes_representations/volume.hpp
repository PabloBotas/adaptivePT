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

    Patient_Volume_t(const std::string& file,
                     const Patient_Volume_t::Source_type& source);
    Patient_Volume_t(const std::string& file,
                     const Patient_Volume_t::Source_type& source,
                     const unsigned int& nx, const unsigned int& ny, const unsigned int& nz,
                     const float& dx, const float& dy, const float& dz);
    Patient_Volume_t(const float* src,
                     const CT_Dims_t& dims);
    Patient_Volume_t(const CT_Dims_t& dims);

    void setDimsFromPatient(const CT_Dims_t& pat_ct);
    void setVoxels(unsigned int x, unsigned int y, unsigned int z);
    void setSpacing(float x, float y, float z);
    void output(std::string outfile, std::string out_type);

    Patient_Volume_t::Source_type source_type;
    std::string file;
    Vector_t<unsigned int> n;
    Vector_t<float> d;
    Vector_t<float> origin;
    Vector_t<float> imgCenter;
    unsigned int nElements;

    std::vector<float> hu;
private:
    Vector_t<unsigned int> original_n;
    Vector_t<float> original_d;
    Vector_t<float> original_origin;
    Vector_t<float> original_imgCenter;
    void read_volume();
    void import_from_metaimage(const float* data);
    void export_to_metaimage(std::string f);
    void consolidate_originals();
};

#endif
