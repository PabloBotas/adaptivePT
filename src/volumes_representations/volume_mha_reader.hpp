#ifndef __MHA_VOLUMES_HPP__
#define __MHA_VOLUMES_HPP__

#include "special_types.hpp"

#include <string>
#include <vector>
#include <fstream>

class Mha_reader_t
{
public:
    Mha_reader_t();
    Mha_reader_t(std::string file);

    std::string file;
    Vector3_t<unsigned int> dim;
    unsigned int nElements;
    Vector3_t<double> spacing;
    Vector3_t<double> origin;
    std::vector<double> transform_matrix;
    unsigned short nb;
    unsigned short type_id;

    std::vector<short> data;


    // Header variables
    std::string ObjectType;
    unsigned int NDims;
    bool BinaryData;
    bool BinaryDataByteOrderMSB;
    bool CompressedData;
    std::vector<int> TransformMatrix;
    std::vector<int> CenterOfRotation;
    std::string AnatomicalOrientation;
    std::string ElementType;
    std::string ElementDataFile;
    unsigned int ElementNumberOfChannels;

private:
    void read_file();
    void read_header();
    void read_body();
    template<class T>
    void read_body_to_vector(std::vector<T>& temp, std::ifstream& stream);
    template<class T>
    T getHeaderValue(std::ifstream &stream, std::string key, T default_value = T());
    template<class T>
    std::vector<T> getHeaderVector(std::ifstream &stream, std::string key, const int n);

    unsigned int default_ElementNumberOfChannels = 1;

};

#endif