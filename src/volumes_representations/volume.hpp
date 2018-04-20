#ifndef __PATIENT_VOLUME_HPP__
#define __PATIENT_VOLUME_HPP__

#include "special_types.hpp"

#include <string>
#include <vector>
#include <fstream>


struct Volume_metadata_t
{
    Vector3_t<unsigned int> n;
    Vector3_t<float> d;
    Vector3_t<float> origin;
    Vector3_t<float> imgCenter;
    unsigned int nElements;
    void print() const {
        std::cout << "Vol metadata:" << std::endl;
        std::cout << "\tn: " << n.x << " " << n.y << " " << n.z << std::endl;
        std::cout << "\td: " << d.x << " " << d.y << " " << d.z << std::endl;
        std::cout << "\torigin: " << origin.x << " " << origin.y << " " << origin.z << std::endl;
        std::cout << "\timgCenter: " << imgCenter.x << " " << imgCenter.y << " " << imgCenter.z << std::endl;
        std::cout << "\tnElements: " << nElements << std::endl;
    }
};

class Volume_t
{
public:
    enum Source_type {
        CTVOLUME,
        DOSE,
        MHA
    };

    Volume_t(const std::string& file, const Volume_t::Source_type& source);
    Volume_t(const Volume_metadata_t& meta, bool long_data_ = false);
    Volume_t(const CT_Dims_t& dims, const bool long_data_ = false);

    void freeMemory();
    void setDims(const Volume_metadata_t& dims);
    void setDims(const CT_Dims_t& pat_ct, const bool interpolate = false);
    void setVoxels(unsigned int x, unsigned int y, unsigned int z);
    void setSpacing(float x, float y, float z);
    void setInternalFlag();
    void setExternalFlag();
    void normalize(float ref);
    void scale(float ref);
    void output(std::string outfile);
    void output(std::string outfile, const CT_Dims_t& dims);
    void mha_to_int_coordinates();
    void int_to_mha_coordinates();
    uint count_above_thres(float thres);
    Volume_metadata_t getMetadata() const;

    Volume_t::Source_type source_type;
    std::string file;
    Vector3_t<unsigned int> n;
    Vector3_t<float> d;
    Vector3_t<float> origin;
    Vector3_t<float> imgCenter;
    unsigned int nElements;

    std::vector<float> data;
    std::vector<unsigned long long int> long_data;
private:
    Vector3_t<unsigned int> original_n;
    Vector3_t<float> original_d;
    Vector3_t<float> original_origin;
    Vector3_t<float> original_imgCenter;
    Vector3_t<float> original_mha_origin;
    void read_volume();
    template <class T>
    void import_from_metaimage(const std::vector<T>& data);
    void export_binary(std::string f);
    void export_header_metaimage(std::string f, std::string ref_file = "LOCAL");
    void export_binary_metaimage(std::string f, std::ios::openmode mode = std::ios::out);
    void consolidate_originals();
    void interpolate_to_geometry(Volume_t& pat, const CT_Dims_t& pat_ct, const float extrapolationValue = 0); 
    void interpolate_to_geometry(const CT_Dims_t& pat_ct, const float extrapolationValue = 0); 
    // void interpolate_to_geometry(const Volume_metadata_t& pat_ct, const float extrapolationValue = 0); 
    void do_interpolate(std::vector<float>& dest, const CT_Dims_t& pat_ct, const float extrapolationValue = 0); 
    bool high_precision = false;
    bool internal_coords;
};

#endif
