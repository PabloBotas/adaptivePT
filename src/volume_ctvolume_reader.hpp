#ifndef __VOLUME_CTVOLUME_READER_HPP__
#define __VOLUME_CTVOLUME_READER_HPP__

#include <string>
#include <vector>
#include <fstream>


class Ctvolume_reader_t
{
public:
    Ctvolume_reader_t();
    Ctvolume_reader_t(std::string file);

    unsigned int nElements;
    std::vector<short> hu;

private:
    void read_file(std::string f);
};

#endif