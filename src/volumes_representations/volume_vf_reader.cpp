#include "volume_vf_reader.hpp"
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

Vf_reader_t::Vf_reader_t(std::string f) : file(f)
{
    read_file();
}


void Vf_reader_t::read_file()
{
    read_header();
    read_body();
}


void Vf_reader_t::read_header()
{
    std::ifstream stream(file);
    utils::check_fs(stream, file, std::string());

    // Get header values
    std::string dummy_line;
    std::getline(stream, dummy_line);
    std::getline(stream, dummy_line);
    std::getline(stream, dummy_line);
    std::getline(stream, dummy_line);
    std::getline(stream, dummy_line);
    std::getline(stream, dummy_line);
    origin = getHeaderVector<float>(stream, "Offset", 3);
    std::getline(stream, dummy_line);
    std::getline(stream, dummy_line);
    d = getHeaderVector<float>(stream, "ElementSpacing", 3);
    n = getHeaderVector<unsigned int>(stream, "DimSize", 3);
    std::getline(stream, dummy_line);
    std::string ElType = getHeaderValue<std::string>(stream, "ElementType");
    std::string ElementDataFile = getHeaderValue<std::string>(stream, "ElementDataFile");

    nElements = n.x*n.y*n.z;

    if (ElType.compare("MET_FLOAT") == 0) {
        nb = sizeof(float);
        type_id = 3;
    } else {
        std::cerr << "Unrecognized type \"" << ElType;
        std::cerr << "\" in MHA header of file " << file << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "MHA file: " << file << std::endl;
    std::cout << "    - dimensions: ";
    std::cout << n.x << ", " << n.y << ", " << n.z << " (" << nElements << ")" << std::endl;
    std::cout << "    - origin:     ";
    std::cout << origin.x << ", " << origin.y << ", " << origin.z << std::endl;
    std::cout << "    - spacing:    ";
    std::cout << d.x << ", " << d.y << ", " << d.z << std::endl;
    std::cout << "    - data type:  ";
    std::cout << ElType << " (" << nb << " bytes)" << std::endl;
}


void Vf_reader_t::read_body()
{
    unsigned int bytes_to_read = nElements*nb*3;
    std::ifstream stream(file, std::ios::binary | std::ios::in | std::ios::ate);
    unsigned int end = stream.tellg();
    unsigned int header_size = end - bytes_to_read;
    stream.seekg(header_size);

    data.resize(nElements);
    switch (type_id)
    {
        case 3:
        {
            // std::cout << "INPUT IMAGE FORMAT IS FLOAT" << std::endl;
            std::vector<float> temp_x(nElements);
            std::vector<float> temp_y(nElements);
            std::vector<float> temp_z(nElements);
            for (unsigned int i = 0; i < nElements; ++i) {
                float temp_x;
                float temp_y;
                float temp_z;
                stream.read(reinterpret_cast<char*>(&temp_x), nb);
                // stream.read(reinterpret_cast<char*>(&temp_z), nb);
                stream.read(reinterpret_cast<char*>(&temp_y), nb);
                stream.read(reinterpret_cast<char*>(&temp_z), nb);
                // stream.read(reinterpret_cast<char*>(&temp_x), nb);
                // data.at(i).x = temp_x[i]/10;
                // data.at(i).y = -temp_y[i]/10;
                // data.at(i).z = -temp_z[i]/10;
                data.at(i).x = temp_x;
                data.at(i).y = temp_y;
                data.at(i).z = temp_z;

                // if (temp_x != 0 || temp_y != 0 || temp_z != 0) {
                //     std::cout << temp_x << " " << temp_y << " " << temp_z << std::endl;
                // }
            }
            break;
        }
        default:
            break;
    }

    std::cout << "Bytes read: " << bytes_to_read << std::endl;
}

void Vf_reader_t::to_int_coordinates()
{
    std::swap(d.x, d.z);
    std::swap(n.x, n.z);
    std::swap(origin.x, origin.z);
    origin.x *= -1;

    Array3<float> temp = data;
    // std::copy(data.begin(), data.end(), temp.begin());
    for (size_t i=0; i<nElements; ++i) {
        temp.at(i) /= 10;
        std::swap(temp.at(i).x, temp.at(i).z);
        temp.at(i).x *= 1;
        temp.at(i).y *= -1;
        temp.at(i).z *= -1;

        int ix = i/(n.y*n.z) % n.x;
        int iy = i/n.z % n.y;
        int iz = i % n.z;

        ix = n.x-ix-1;

        int idx = iz + iy*n.z + ix*n.z*n.y;
        data.at(idx) = temp.at(i);
    }

    origin.x /= 10;
    origin.y /= 10;
    origin.z /= 10;
    d.x /= 10;
    d.y /= 10;
    d.z /= 10;

    // std::ofstream fout("vf_internal.dat", std::ios::out | std::ios::binary);
    // for (uint i = 0; i < data.size(); ++i) {
    //     float d = std::sqrt(data.at(i).x*data.at(i).x + data.at(i).y*data.at(i).y + data.at(i).z*data.at(i).z);
    //     fout.write((char*)&d, sizeof(float));
    // }
}


template<class T>
T Vf_reader_t::getHeaderValue(std::ifstream &stream, std::string key, T default_value)
{
    size_t place = stream.tellg();

    std::string line;
    std::getline(stream, line);
    std::string field, dummy;
    T out;
    if (typeid(T) == typeid(bool)) {
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);
        std::istringstream ss(line);
        ss >> field >> dummy >> std::boolalpha >> out;
    } else {
        std::istringstream ss(line);
        ss >> field >> dummy >> out;
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
std::vector<T> Vf_reader_t::getHeaderVector(std::ifstream &stream, std::string key, const int n)
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
