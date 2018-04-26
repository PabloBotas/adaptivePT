#include "volume_mha_reader.hpp"
#include "special_types.hpp"
#include "utils.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <limits>
#include <typeinfo>

Mha_reader_t::Mha_reader_t(std::string f) : file(f)
{
    read_file();
}


void Mha_reader_t::read_file()
{
    std::ifstream stream(file);
    utils::check_fs(stream, file, std::string());
    read_header(stream);
    read_body(stream);
}


void Mha_reader_t::read_header(std::ifstream& stream)
{
    // Get header values
    ObjectType              = getHeaderValue<std::string>(stream, "ObjectType");
    NDims                   = getHeaderValue<unsigned int>(stream, "NDims");
    BinaryData              = getHeaderValue<bool>(stream, "BinaryData");
    BinaryDataByteOrderMSB  = getHeaderValue<bool>(stream, "BinaryDataByteOrderMSB");
    CompressedData          = getHeaderValue<bool>(stream, "CompressedData");
    transform_matrix        = getHeaderVector<float>(stream, "TransformMatrix", 9);
    origin                  = getHeaderVector<float>(stream, "Offset", 3);
    CenterOfRotation        = getHeaderVector<int>(stream, "CenterOfRotation", 9);
    AnatomicalOrientation   = getHeaderValue<std::string>(stream, "AnatomicalOrientation");
    spacing                 = getHeaderVector<float>(stream, "ElementSpacing", 3);
    std::string ITK_Input   = getHeaderValue<std::string>(stream, "ITK_InputFilterName", "");
    dim                     = getHeaderVector<unsigned int>(stream, "DimSize", 3);
    ElementNumberOfChannels = getHeaderValue<unsigned int>(stream, "ElementNumberOfChannels",
                                                           default_ElementNumberOfChannels);
    std::string ElType      = getHeaderValue<std::string>(stream, "ElementType");
    ElementDataFile         = getHeaderValue<std::string>(stream, "ElementDataFile");

    origin.x /= 10;
    origin.y /= 10;
    origin.z /= 10;
    spacing.x /= 10;
    spacing.y /= 10;
    spacing.z /= 10;
    nElements = dim.x*dim.y*dim.z;

    if (ElType.compare("MET_SHORT") == 0) {
        nb = sizeof(short);
        type_id = 1;
    } else if (ElType.compare("MET_USHORT") == 0) {
        nb = sizeof(unsigned short);
        type_id = 2;
    } else if (ElType.compare("MET_FLOAT") == 0) {
        nb = sizeof(float);
        type_id = 3;
    } else if (ElType.compare("MET_UCHAR") == 0) {
        nb = sizeof(unsigned char);
        type_id = 4;
    } else {
        std::cerr << "Unrecognized type \"" << ElType;
        std::cerr << "\" in MHA header of file " << file << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string str("MHA file: ");
    str += file + "\n";
    str += "    - dimensions: " + std::to_string(dim.x) + ", ";
    str += std::to_string(dim.y) + ", " + std::to_string(dim.z);
    str += " (" + std::to_string(nElements) + ")\n";

    str += "    - origin:     " + std::to_string(origin.x) + ", ";
    str += std::to_string(origin.y) + ", " + std::to_string(origin.z) + "\n";

    str += "    - spacing:    " + std::to_string(spacing.x) + ", ";
    str += std::to_string(spacing.y) + ", " + std::to_string(spacing.z) + "\n";

    str += "    - data type:  " + ElType + " (" + std::to_string(nb) + ")\n";

    std::cout << str;
}


void Mha_reader_t::read_body(std::ifstream& stream)
{
    data.resize(nElements);
    switch (type_id) {
        case 1:
        {
            // std::cout << "INPUT IMAGE FORMAT IS SHORT" << std::endl;
            std::vector<short> temp(nElements);
            read_body_to_vector(data, stream);
            break;
        }
        case 2:
        {
            // std::cout << "INPUT IMAGE FORMAT IS UNSIGNED SHORT" << std::endl;
            std::vector<unsigned short> temp(nElements);
            read_body_to_vector(temp, stream);
            data.assign(temp.begin(), temp.end());
            break;
        }
        case 3:
        {
            // std::cout << "INPUT IMAGE FORMAT IS FLOAT" << std::endl;
            std::vector<float> temp(nElements);
            read_body_to_vector(temp, stream);
            data.assign(temp.begin(), temp.end());
            break;
        }
        case 4:
        {
            // std::cout << "INPUT IMAGE FORMAT IS UNSIGNED CHAR" << std::endl;
            std::vector<unsigned char> temp(nElements);
            read_body_to_vector<unsigned char>(temp, stream);
            data.assign(temp.begin(), temp.end());
            break;
        }
        default:
        {
            std::cerr << "ERROR! Unrecognized MHA data type of file " << file << std::cerr;
            exit(EXIT_FAILURE);
        }
    }
}

template<class T>
void Mha_reader_t::read_body_to_vector(std::vector<T>& temp, std::ifstream& stream)
{
    if (ElementNumberOfChannels > 1) {
        for (unsigned int i = 0; i < nElements; ++i) {
            stream.read(reinterpret_cast<char*>(&temp[i]), nb);
            stream.ignore((ElementNumberOfChannels-1)*nb);
        }
    } else {
        stream.read(reinterpret_cast<char*>(&temp[0]), nElements*nb);
    }
}


template<class T>
T Mha_reader_t::getHeaderValue(std::ifstream &stream, std::string key, T default_value)
{
    size_t place = stream.tellg();

    std::string line;
    std::getline(stream, line);
    std::string field, equal_sign;
    T out;
    if (typeid(T) == typeid(bool)) {
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);
        std::istringstream ss(line);
        ss >> field >> equal_sign >> std::boolalpha >> out;
    } else {
        std::istringstream ss(line);
        ss >> field >> equal_sign >> out;
    }

    std::transform(key.begin(), key.end(), key.begin(), ::tolower);
    std::transform(field.begin(), field.end(), field.begin(), ::tolower);
    if (field != key) {
        stream.seekg(place);
        return default_value;
    }

    return out;
}

template<class T>
std::vector<T> Mha_reader_t::getHeaderVector(std::ifstream &stream, std::string key, const int n)
{
    size_t place = stream.tellg();
    std::string line;
    std::getline(stream, line);

    std::istringstream ss(line);
    std::string field, dummy;
    ss >> field >> dummy;

    std::transform(key.begin(), key.end(), key.begin(), ::tolower);
    std::transform(field.begin(), field.end(), field.begin(), ::tolower);
    if (field != key) {
        stream.seekg(place);
        return std::vector<T>();
    }

    std::vector<T> out(n);
    for (int i = 0; i < n; i++) {
        ss >> out[i];
    }

    return out;
}
