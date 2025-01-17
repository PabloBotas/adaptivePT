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
    Vector3_t<unsigned int> n;
    unsigned int nElements;
    Vector3_t<float> d;
    Vector3_t<float> origin;
    unsigned short nb;
    unsigned short type_id;

    Array3<float> data;

    void to_int_coordinates();

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