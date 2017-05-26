#include "volume.hpp"
#include "volume_mha_reader.hpp"
#include "volume_ctvolume_reader.hpp"
#include "special_types.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

#define CM2MM 10

Patient_Volume_t::Patient_Volume_t(const std::string& f,
                                   const Patient_Volume_t::Source_type& s)
{
    file = f;
    source_type = s;

    read_volume();
    imgCenter.x = (origin.x - d.x/2) + 1/2*d.x*n.x;
    imgCenter.y = (origin.y - d.y/2) + 1/2*d.y*n.y;
    imgCenter.z = (origin.z - d.z/2) + 1/2*d.z*n.z;
    consolidate_originals();
}

Patient_Volume_t::Patient_Volume_t(const std::string& f,
                                   const Patient_Volume_t::Source_type& s,
                                   const unsigned int& nx,
                                   const unsigned int& ny,
                                   const unsigned int& nz,
                                   const float& dx, const float& dy, const float& dz)
{
    file = f;
    source_type = s;

    read_volume();
    setVoxels(nx, ny, nz);
    setSpacing(dx, dy, dz);
    consolidate_originals();
}

Patient_Volume_t::Patient_Volume_t(const CT_Dims_t& dims)
{
    nElements = dims.total;
    hu.resize(nElements);
    setDimsFromPatient(dims);
    consolidate_originals();
}

Patient_Volume_t::Patient_Volume_t(const float* src,
                                   const CT_Dims_t& dims) :
                                   Patient_Volume_t(dims)
{
    std::copy( src, src+nElements, hu.begin() );
}

void Patient_Volume_t::consolidate_originals()
{
    original_n = n;
    original_d = d;
    original_origin = origin;
    original_imgCenter = imgCenter;
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
            import_from_metaimage((float*)reader.data.data());
            break;
        }
    }
}

void Patient_Volume_t::import_from_metaimage(const float* data)
{
    for (size_t k = 0; k < n.z; k++) {
        for (size_t j = 0; j < n.y; j++) {
            for (size_t i = 0; i < n.x; i++) {
                size_t in  = i + j*n.x +         k*n.x*n.y;
                size_t out = i + j*n.x + (n.z-k-1)*n.x*n.y;
                hu[out] = (float)data[in];
            }
        }
    }
}

void Patient_Volume_t::export_to_metaimage(std::string f)
{
    std::ofstream ofs;
    ofs.open (f, std::ios::out | std::ios::binary);
    if( !ofs.is_open() )
    {
        std::cerr << "Can not open file " << f << " to write results." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::vector<float> temp(nElements);
    for (size_t k = 0; k < n.x; k++) {
        for (size_t j = 0; j < n.y; j++) {
            for (size_t i = 0; i < n.z; i++) {
                size_t in  = i + j*n.z +         k*n.z*n.y;
                size_t out = i + j*n.z + (n.x-k-1)*n.z*n.y;
                temp[out] = hu[in];
            }
        }
    }
    ofs.write (reinterpret_cast<char*>(temp.data()), nElements*sizeof(float));
    ofs.close();
}

void Patient_Volume_t::output(std::string outfile, std::string out_type)
{
    if (out_type != "mhd")
    {
        std::cerr << "Output type not supported, supported type: mhd" << std::endl;
    }

    std::string file_basename = outfile.substr(0, outfile.find_last_of("/"));
    std::string file_noext = outfile.substr(0, outfile.find_last_of("."));
    std::string mhd_header_file = file_noext + ".mhd";

    // Output header file
    std::ofstream header(mhd_header_file, std::ios::out);
    if( !header.is_open() )
    {
        std::cerr << "Can not open file " << mhd_header_file << " to write results." << std::endl;
        exit(EXIT_FAILURE);
    }
    header << "ObjectType = Image\n";
    header << "NDims = 3\n";
    header << "BinaryData = True\n";
    header << "BinaryDataByteOrderMSB = False\n";
    header << "TransformMatrix = 1 0 0 0 1 0 0 0 1\n";
    header << "Offset = " << CM2MM*origin.z << " " << CM2MM*origin.y << " " << CM2MM*origin.x << "\n";
    header << "CenterOfRotation = 0 0 0\n";
    header << "ElementSpacing = " << CM2MM*d.z << " " << CM2MM*d.y << " " << CM2MM*d.x << "\n";
    header << "DimSize = " << n.z << " " << n.y << " " << n.x << "\n";
    header << "AnatomicalOrientation = RAI\n";
    header << "ElementType = MET_FLOAT\n";
    header << "ElementDataFile = " << file_basename << "\n";
    header.close();

    std::cout << "Offset = " << CM2MM*origin.z << " " << CM2MM*origin.y << " " << CM2MM*origin.x << "\n";
    std::cout << "ElementSpacing = " << CM2MM*d.z << " " << CM2MM*d.y << " " << CM2MM*d.x << "\n";
    std::cout << "DimSize = " << n.z << " " << n.y << " " << n.x << "\n";

    // Output binary
    export_to_metaimage(outfile);
}

void Patient_Volume_t::setDimsFromPatient(const CT_Dims_t& pat_ct)
{
    setVoxels(pat_ct.n.x, pat_ct.n.y, pat_ct.n.z);
    setSpacing(pat_ct.d.x, pat_ct.d.y, pat_ct.d.z);
    origin.x = pat_ct.offset.x + 0.5*pat_ct.d.x + pat_ct.isocenter.x;
    origin.y = pat_ct.offset.y + 0.5*pat_ct.d.y + pat_ct.isocenter.y;
    origin.z = pat_ct.offset.z + 0.5*pat_ct.d.z + pat_ct.isocenter.z;
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
