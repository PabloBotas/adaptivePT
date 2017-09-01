#ifndef __VOLUME_DOSE_READER_HPP__
#define __VOLUME_DOSE_READER_HPP__

#include <string>
#include <vector>
#include <fstream>


class Dose_reader_t
{
public:
    Dose_reader_t();
    Dose_reader_t(std::string file);

    unsigned int nElements;
    std::vector<float> data;

private:
    void read_file(std::string f);
};

#endif