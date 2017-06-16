#ifndef __PATIENT_VOLUME_HPP__
#define __PATIENT_VOLUME_HPP__

#include "special_types.hpp"

#include <string>
#include <vector>
#include <fstream>

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

    void setDims(const CT_Dims_t& pat_ct, const bool interpolate = false);
    void setVoxels(unsigned int x, unsigned int y, unsigned int z);
    void setSpacing(float x, float y, float z);
    void output(std::string outfile, std::string out_type);
    void output(std::string outfile, std::string out_type, const CT_Dims_t& dims);

    Patient_Volume_t::Source_type source_type;
    std::string file;
    Vector3_t<unsigned int> n;
    Vector3_t<float> d;
    Vector3_t<float> origin;
    Vector3_t<float> imgCenter;
    unsigned int nElements;

    std::vector<float> hu;
private:
    Vector3_t<unsigned int> original_n;
    Vector3_t<float> original_d;
    Vector3_t<float> original_origin;
    Vector3_t<float> original_imgCenter;
    void read_volume();
    template <class T>
    void import_from_metaimage(const std::vector<T>& data);
    void export_binary(std::string f);
    void export_header_metaimage(std::string f, std::string ref_file = "LOCAL");
    void export_binary_metaimage(std::string f, std::ios::openmode mode = std::ios::out);
    void consolidate_originals();
    void interpolate_to_geometry(Patient_Volume_t& pat, const CT_Dims_t& pat_ct, const float extrapolationValue = 0);
    void interpolate_to_geometry(const CT_Dims_t& pat_ct, const float extrapolationValue = 0);
    void do_interpolate(std::vector<float>& dest, const CT_Dims_t& pat_ct, const float extrapolationValue = 0);

};

#endif
