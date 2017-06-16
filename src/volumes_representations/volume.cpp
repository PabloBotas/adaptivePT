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
    setDims(dims);
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
            import_from_metaimage<short>(reader.data);
            break;
        }
    }
}

template <class T>
void Patient_Volume_t::import_from_metaimage(const std::vector<T>& vec)
{
    if (hu.size() != vec.size())
    {
        std::cerr << "ERROR! in file " << __FILE__ << ", line " << __LINE__;
        std::cerr << "\nSizes are incompatible: " << hu.size() << " != " << vec.size();
        std::cerr << std::endl;
        exit(EXIT_FAILURE);
    }
    for (size_t k = 0; k < n.z; k++) {
        for (size_t j = 0; j < n.y; j++) {
            for (size_t i = 0; i < n.x; i++) {
                size_t in  = i + j*n.x +         k*n.x*n.y;
                size_t out = i + j*n.x + (n.z-k-1)*n.x*n.y;
                hu.at(out) = (float)vec.at(in);
            }
        }
    }
}
template void Patient_Volume_t::import_from_metaimage<short>(const std::vector<short>& data);

void Patient_Volume_t::export_binary_metaimage(std::string f,
                                               std::ios::openmode other)
{
    std::ios::openmode mode = std::ios::out | std::ios::binary | other;

    std::ofstream ofs;
    ofs.open (f, mode);
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

void Patient_Volume_t::export_binary(std::string f)
{
    std::ofstream ofs;
    ofs.open (f, std::ios::out | std::ios::binary);
    if( !ofs.is_open() )
    {
        std::cerr << "Can not open file " << f << " to write results." << std::endl;
        exit(EXIT_FAILURE);
    }

    ofs.write (reinterpret_cast<char*>(hu.data()), nElements*sizeof(float));
    ofs.close();
}

void Patient_Volume_t::output(std::string outfile, std::string out_type)
{
    if (out_type == "mhd")
    {
        std::string file_basename = outfile.substr(0, outfile.find_last_of("/"));
        std::string file_noext = outfile.substr(0, outfile.find_last_of("."));
        std::string mhd_header_file = file_noext + ".mhd";
        export_header_metaimage(mhd_header_file, file_basename);
        export_binary_metaimage(outfile);
    }
    else if(out_type == "mha")
    {
        export_header_metaimage(outfile);
        export_binary_metaimage(outfile, std::ios::app);
    }
    else if(out_type == "bin" || out_type == "binary")
    {
        export_binary(outfile);
    }
    else
        std::cerr << "Output type not supported, supported type: mhd" << std::endl;

    // Export file
}

void Patient_Volume_t::export_header_metaimage(std::string outfile, std::string ref_file)
{
    // Output header file
    std::ofstream header(outfile, std::ios::out);
    if( !header.is_open() )
    {
        std::cerr << "Can not open file " << outfile << " to write results." << std::endl;
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
    header << "ElementDataFile = " << ref_file << "\n";
    header.close();

    std::cout << "Offset = " << CM2MM*origin.z << " " << CM2MM*origin.y << " " << CM2MM*origin.x << "\n";
    std::cout << "ElementSpacing = " << CM2MM*d.z << " " << CM2MM*d.y << " " << CM2MM*d.x << "\n";
    std::cout << "DimSize = " << n.z << " " << n.y << " " << n.x << "\n";
}

void Patient_Volume_t::output(std::string outfile, std::string out_type, const CT_Dims_t& dims)
{
    Patient_Volume_t temp(dims);
    this->interpolate_to_geometry(temp, dims, 0);
    temp.output(outfile, out_type);
}

void Patient_Volume_t::setDims(const CT_Dims_t& pat_ct, const bool interpolate)
{
    setVoxels(pat_ct.n.x, pat_ct.n.y, pat_ct.n.z);
    setSpacing(pat_ct.d.x, pat_ct.d.y, pat_ct.d.z);
    origin.x = pat_ct.offset.x + 0.5*pat_ct.d.x + pat_ct.isocenter.x;
    origin.y = pat_ct.offset.y + 0.5*pat_ct.d.y + pat_ct.isocenter.y;
    origin.z = pat_ct.offset.z + 0.5*pat_ct.d.z + pat_ct.isocenter.z;

    if (interpolate)
        interpolate_to_geometry(pat_ct);
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

void Patient_Volume_t::interpolate_to_geometry(const CT_Dims_t& pat_ct,
                                               const float extrapolationValue)
{
    do_interpolate(hu, pat_ct, extrapolationValue);
}

void Patient_Volume_t::interpolate_to_geometry(Patient_Volume_t& pat,
                                               const CT_Dims_t& pat_ct,
                                               const float extrapolationValue)
{
    do_interpolate(pat.hu, pat_ct, extrapolationValue);
}

void Patient_Volume_t::do_interpolate(std::vector<float>& dest,
                                      const CT_Dims_t& pat_ct,
                                      const float extrapolationValue)
{
    for (size_t tid = 0; tid < pat_ct.total; tid++)
    {
        size_t zInd = tid%pat_ct.n.z;
        size_t yInd = (tid/pat_ct.n.z)%pat_ct.n.y;
        size_t xInd = tid/pat_ct.n.z/pat_ct.n.y;

        // 3d volumetric interpolation
        float fraction, interpolatedValue = 0.0f, mass = 0.0f;
        Vector3_t<float> inf, sup;
        //calculate the position of the scoring grid boundary
        inf.x = xInd*pat_ct.d.x + pat_ct.offset.x;
        inf.y = yInd*pat_ct.d.y + pat_ct.offset.y;
        inf.z = zInd*pat_ct.d.z + pat_ct.offset.z;
        sup.x = xInd*(pat_ct.d.x + 1) + pat_ct.offset.x;
        sup.y = yInd*(pat_ct.d.y + 1) + pat_ct.offset.y;
        sup.z = zInd*(pat_ct.d.z + 1) + pat_ct.offset.z;

        //calculate the nearest CT grid walls
        Vector3_t<int> infCT, supCT;
        infCT.x = (int) floorf(inf.x/this->d.x);
        infCT.y = (int) floorf(inf.y/this->d.y);
        infCT.z = (int) floorf(inf.z/this->d.z);
        supCT.x = (int) floorf(sup.x/this->d.x);
        supCT.y = (int) floorf(sup.y/this->d.y);
        supCT.z = (int) floorf(sup.z/this->d.z);

        //loop over the CT grid voxels which overlap with the current interpolated grid voxel
        for(int zvox = infCT.z; zvox <= supCT.z; zvox++) {
            for(int yvox = infCT.y; yvox <= supCT.y; yvox++) {
                for(int xvox = infCT.x; xvox <= supCT.x; xvox++) {

                    float ctValue;

                    if (! (xvox < 0 || xvox >= (int) this->n.x
                        || yvox < 0 || yvox >= (int) this->n.y
                        || zvox < 0 || zvox >= (int) this->n.z)) {
                        // Inside CT grid
                        ctValue = hu.at(zvox + yvox*this->n.z + xvox*this->n.z*this->n.y);
                    } else {
                        // Outside CT grid
                        ctValue = extrapolationValue;
                    }

                    Vector3_t<float> fractions;
                    // determine the overlap along x direction
                    if(xvox == infCT.x) {
                        if(xvox == supCT.x)
                            fractions.x = pat_ct.d.x / this->d.x;
                        else
                            fractions.x = ((xvox+1)*this->d.x - inf.x) / this->d.x;
                    } else {
                        if(xvox == supCT.x)
                            fractions.x = (sup.x - xvox*this->d.x) / this->d.x;
                        else
                            fractions.x = 1.0f;
                    }
                    // determine the overlap along y direction
                    if(yvox == infCT.y) {
                        if(yvox == supCT.y)
                            fractions.y = pat_ct.d.y / this->d.y;
                        else
                            fractions.y = ((yvox+1)*this->d.y - inf.y) / this->d.y;
                    } else {
                        if(yvox == supCT.y)
                            fractions.y = (sup.y - yvox*this->d.y) / this->d.y;
                        else
                            fractions.y = 1.0f;
                    }
                    // determine the overlap along z direction
                    if(zvox == infCT.z) {
                        if(zvox == supCT.z)
                            fractions.z = pat_ct.d.z / this->d.z;
                        else
                            fractions.z = ((zvox+1)*this->d.z - inf.z) / this->d.z;
                    } else {
                        if(zvox == supCT.z)
                            fractions.z = (sup.z - zvox*this->d.z) / this->d.z;
                        else
                            fractions.z = 1.0f;
                    }

                    // overlap in 3d
                    fraction = fractions.x * fractions.y * fractions.z;
                    mass += fraction*this->d.x*this->d.y*this->d.z;
                    interpolatedValue += fraction*ctValue*this->d.x*this->d.y*this->d.z;
                }
            }
        }
        if(mass > 0.0) {
            dest.at(tid) = interpolatedValue/mass;
        } else {
            dest.at(tid) = extrapolationValue;
        }
    }
}


