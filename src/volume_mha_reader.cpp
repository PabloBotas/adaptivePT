#include "volume_mha_reader.hpp"
#include "special_types.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <limits>

Mha_reader_t::Mha_reader_t(std::string f) : file(f)
{
    read_file();
}

void Mha_reader_t::read_file()
{
    read_header();
    read_body();
}


void Mha_reader_t::read_header()
{
    std::ifstream stream(file);
    if (!stream.is_open()) {
        std::cerr << "Can't open file: " << file << std::endl;
        return;
    }

    // Get header values
    ObjectType             = getHeaderValue<std::string>(stream);
    NDims                  = getHeaderValue<unsigned int>(stream);
    BinaryData             = getHeaderValue<bool>(stream);
    BinaryDataByteOrderMSB = getHeaderValue<bool>(stream);
    CompressedData         = getHeaderValue<bool>(stream);
    std::vector<float> TransformMatrix = getHeaderVector<float>(stream, 9);
    std::vector<float> Offset          = getHeaderVector<float>(stream, 3);
    CenterOfRotation       = getHeaderVector<int>(stream, 9);
    AnatomicalOrientation  = getHeaderValue<std::string>(stream);
    std::vector<float> ElementSpacing = getHeaderVector<float>(stream, 3);
    std::vector<unsigned int> DimSize = getHeaderVector<unsigned int>(stream, 3);
    std::string ElementType = getHeaderValue<std::string>(stream);
    ElementDataFile         = getHeaderValue<std::string>(stream);

    transform_matrix = TransformMatrix;
    origin           = Offset;
    dim              = DimSize;
    nElements        = dim.x*dim.y*dim.z;
    spacing          = ElementSpacing;

    if (ElementType.compare("MET_SHORT") == 0) {
        nb = 2;
        type_id = 1;
    } else if (ElementType.compare("MET_USHORT") == 0) {
        nb = 2;
        type_id = 2;
    } else if (ElementType.compare("MET_FLOAT") == 0) {
        nb = 4;
        type_id = 3;
    } else {
        std::cerr << "Unrecognized type" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "MHA file: " << file << std::endl;
    std::cout << "    - dimensions: ";
    std::cout << dim.x << ", " << dim.y << ", " << dim.z << " (" << nElements << ")" << std::endl;
    std::cout << "    - origin:     ";
    std::cout << origin.x << ", " << origin.y << ", " << origin.z << std::endl;
    std::cout << "    - spacing:    ";
    std::cout << spacing.x << ", " << spacing.y << ", " << spacing.z << std::endl;
    std::cout << "    - data type:  ";
    std::cout << ElementType << " (" << nb << " bytes)" << std::endl;
}


void Mha_reader_t::read_body()
{
    unsigned int bytes_to_read = nElements*nb;
    std::ifstream stream(file, std::ios::binary | std::ios::in | std::ios::ate);
    unsigned int end = stream.tellg();
    unsigned int header_size = end - bytes_to_read;
    stream.seekg(header_size);

    data.resize(nElements);
    switch (type_id)
    {
        case 1:
        {
            // std::cout << "INPUT IMAGE FORMAT IS SHORT" << std::endl;
            std::vector<short> temp(nElements);
            stream.read(reinterpret_cast<char*>(&data[0]), bytes_to_read);
            break;
        }
        case 2:
        {
            // std::cout << "INPUT IMAGE FORMAT IS UNSIGNED SHORT" << std::endl;
            std::vector<unsigned short> temp(nElements);
            stream.read(reinterpret_cast<char*>(&temp[0]), bytes_to_read);
            data.assign(temp.begin(), temp.end());
            break;
        }
        case 3:
        {
            // std::cout << "INPUT IMAGE FORMAT IS FLOAT" << std::endl;
            std::vector<float> temp(nElements);
            stream.read(reinterpret_cast<char*>(&temp[0]), bytes_to_read);
            data.assign(temp.begin(), temp.end());
            break;
        }
        default:
            break;
    }

    std::cout << "Bytes read: ";
    std::cout << bytes_to_read << std::endl;
}


template<class T>
T Mha_reader_t::getHeaderValue(std::ifstream &stream)
{
    std::string line;
    std::getline(stream, line);
    std::string dummy;
    T out;
    if (typeid(T) == typeid(bool))
    {
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);
        std::istringstream ss(line);
        ss >> dummy >> dummy >> std::boolalpha >> out;
    }
    else
    {
        std::istringstream ss(line);
        ss >> dummy >> dummy >> out;
    }
    return out;
}


template<class T>
std::vector<T> Mha_reader_t::getHeaderVector(std::ifstream &stream, const int n)
{
    std::string line;
    std::getline(stream, line);

    std::istringstream ss(line);
    std::string dummy;
    std::vector<T> out(n);
    ss >> dummy >> dummy;
    for (int i = 0; i < n; i++)
    {
        ss >> out[i];
    }
    return out;
}