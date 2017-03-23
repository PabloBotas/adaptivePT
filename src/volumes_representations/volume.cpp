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
    imgCenter.x = (origin.x - d.x/2) + 1/2*d.x*n.x;
    imgCenter.y = (origin.y - d.y/2) + 1/2*d.y*n.y;
    imgCenter.z = (origin.z - d.z/2) + 1/2*d.z*n.z;
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
            break;
        }
        case Source_type::MHA:
        {
            Mha_reader_t reader(file);

            nElements = reader.nElements;
            n = reader.dim;
            d = reader.spacing;
            origin = reader.origin;
            hu.resize(nElements);
            for (size_t k = 0; k < reader.dim.z; k++)
            {
                for (size_t j = 0; j < reader.dim.y; j++)
                {
                    for (size_t i = 0; i < reader.dim.x; i++)
                    {
                        size_t in  = i + j*reader.dim.x +                  k*reader.dim.x*reader.dim.y;
                        size_t out = i + j*reader.dim.x + (reader.dim.z-k-1)*reader.dim.x*reader.dim.y;
                        hu[out] = (float)reader.data[in];
                    }
                }
            }
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

