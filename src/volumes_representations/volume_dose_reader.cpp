#include "volume_dose_reader.hpp"

#include "utils.hpp"

#include <iostream>
#include <string>
#include <vector>


Dose_reader_t::Dose_reader_t(std::string f)
{
    read_file(f);
}

void Dose_reader_t::read_file(std::string file)
//    Reads voxel geometry from an input file
{

    std::cout << "Reading " << file << std::endl;

    //  open file
    std::ifstream stream(file, std::ios::binary | std::ios::ate);
    utils::check_fs(stream, file, "to read planning CT.");
    nElements = stream.tellg()/sizeof(float);
    size_t bytes_to_read = nElements*sizeof(float);
    stream.seekg(0, std::ios::beg);
    data.resize(nElements);
    stream.read(reinterpret_cast<char*>(&data[0]), bytes_to_read);
}
