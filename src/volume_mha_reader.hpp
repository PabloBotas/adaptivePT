#ifndef __MHA_VOLUMES_HPP__
#define __MHA_VOLUMES_HPP__

#include <string>
#include <vector>
#include <fstream>

template<class T>
struct Vector_t
{
    Vector_t()
    {
        x = 0;
        y = 0;
        z = 0;
    }

    Vector_t(const std::vector<T> &v)
    {
        x = v.at(0);
        y = v.at(1);
        z = v.at(2);
    }

    Vector_t(const Vector_t &obj)
    {
        x = obj.x;
        y = obj.y;
        z = obj.z;
    }

    T x;
    T y;
    T z;
};


class Mha_t
{
public:
    Mha_t();
    Mha_t(std::string file);

    std::string file;
    Vector_t<unsigned int> dim;
    unsigned int nElements;
    Vector_t<double> spacing;
    Vector_t<double> origin;
    std::vector<double> transform_matrix;
    unsigned short nb;
    unsigned short type_id;

    std::vector<short> data;


    // Header variables
    std::string ObjectType;
    unsigned int NDims;
    bool BinaryData;
    bool BinaryDataByteOrderMSB;
    bool CompressedData;
    std::vector<int> TransformMatrix;
    std::vector<int> CenterOfRotation;
    std::string AnatomicalOrientation;
    std::string ElementType;
    std::string ElementDataFile;

private:
    void read_file();
    void read_header();
    void read_body();
    template<class T>
    T getHeaderValue(std::ifstream &stream);
    template<class T>
    std::vector<T> getHeaderVector(std::ifstream &stream, const int n);

};

#endif