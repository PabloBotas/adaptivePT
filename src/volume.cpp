#include "volume.hpp"
#include "volume_mha_reader.hpp"
#include "volume_ctvolume_reader.hpp"
#include "special_types.hpp"

#include <string>
#include <vector>
#include <algorithm>

Patient_Volume_t::Patient_Volume_t()
{
}

Patient_Volume_t::Patient_Volume_t(std::string f,
                                   Patient_Volume_t::Source_type s)
{
    file = f;
    source_type = s;

    read_volume();
}

Patient_Volume_t::Patient_Volume_t(std::string f,
                                   Patient_Volume_t::Source_type s,
                                   unsigned int nx, unsigned int ny, unsigned int nz,
                                   float dx, float dy, float dz)
{
    file = f;
    source_type = s;

    read_volume();
    setVoxels(nx, ny, nz);
    setSpacing(dx, dy, dz);
}

void Patient_Volume_t::read_volume()
{
    switch(source_type)
    {
        case Source_type::CTVOLUME:
        {
            Ctvolume_reader_t reader(file);

            nElements = reader.nElements;
            hu.resize(nElements);
            std::copy( reader.hu.begin(), reader.hu.end(), hu.begin() );
            // hu = new std::vector<float>(reader.hu.begin(), reader.hu.end());
            break;
        }
        case Source_type::MHA:
        {
            Mha_reader_t reader(file);

            nElements = reader.nElements;
            n = reader.dim;
            d = reader.spacing;
            hu.resize(nElements);
            std::copy(reader.data.begin(), reader.data.end(), hu.begin());
            // hu = new std::vector<float>(reader.data.begin(), reader.data.end());
            break;
        }
        case Source_type::MHA_CALIBRATION:
        {
            Mha_reader_t reader(file);
            reader.calibrate(-32678.0, 1);

            nElements = reader.nElements;
            n = reader.dim;
            d = reader.spacing;
            hu.resize(nElements);
            std::copy(reader.data.begin(), reader.data.end(), hu.begin());
            // hu = new std::vector<float>(reader.data.begin(), reader.data.end());
            break;
        }
    }
}

void Patient_Volume_t::setVoxels(unsigned int x, unsigned int y, unsigned int z)
{
    n.x = x;
    n.y = y;
    n.z = z;
}

void Patient_Volume_t::setSpacing(float x, float y, float z)
{
    d.x = x;
    d.y = y;
    d.z = z;
}

