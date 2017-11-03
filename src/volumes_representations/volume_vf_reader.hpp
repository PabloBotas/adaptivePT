#ifndef __VF_VOLUMES_HPP__
#define __VF_VOLUMES_HPP__

#include "special_types.hpp"

#include <string>
#include <vector>
#include <fstream>

class Vf_reader_t
{
public:
    Vf_reader_t();
    Vf_reader_t(std::string file);

    std::string file;
    Vector3_t<unsigned int> dim;
    unsigned int nElements;
    Vector3_t<double> spacing;
    Vector3_t<double> origin;
    unsigned short nb;
    unsigned short type_id;

    Array3<float> data;

private:
    void read_file();
    void read_header();
    void read_body();
    template<class T>
    T getHeaderValue(std::ifstream &stream, std::string key, T default_value = T());
    template<class T>
    std::vector<T> getHeaderVector(std::ifstream &stream, std::string key, const int n);
};

#endif