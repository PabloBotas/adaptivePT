#include "volume_ctvolume_reader.hpp"

#include "utils.hpp"

#include <iostream>
#include <string>
#include <vector>


Ctvolume_reader_t::Ctvolume_reader_t(std::string f)
{
    read_file(f);
}

void Ctvolume_reader_t::read_file(std::string file)
//    Reads voxel geometry from an input file
{

    std::cout << "Reading " << file << std::endl;

    //  open file
    std::ifstream stream(file, std::ios::binary | std::ios::ate);
    utils::check_fs(stream, file, "to read planning CT.");
    nElements = stream.tellg()/sizeof(short);
    size_t bytes_to_read = nElements*sizeof(short);
    stream.seekg(0, std::ios::beg);
    hu.resize(nElements);
    stream.read(reinterpret_cast<char*>(&hu[0]), bytes_to_read);
}
