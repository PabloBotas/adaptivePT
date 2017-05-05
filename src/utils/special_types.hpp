#ifndef __SPECIAL_TYPES_HPP__
#define __SPECIAL_TYPES_HPP__

#include <vector>
#include <string>

template <class T>
class Vector_t
{
public:
    Vector_t();
    Vector_t(const std::vector<T> &v);
    Vector_t(const Vector_t<T> &obj);

    T x;
    T y;
    T z;
};

// Topas_Vector_t ---------------------------------------
template<class T=float, class D=T>
struct Topas_Vector_t
{
    T x;
    T y;
    D z;
};

// BeamAngles_t ---------------------------------------
struct BeamAngles_t
{
    float gantry;
    float couch;
};

// CT_Dims_t ---------------------------------------
struct CT_Dims_t
{
    Topas_Vector_t<> offset;
    Topas_Vector_t<> isocenter;
    Topas_Vector_t<float, std::vector<float> > d;
    Topas_Vector_t<unsigned int, std::vector<unsigned int> > n;
    unsigned int total;
};

// Aperture_Dims_t ---------------------------------------
struct Aperture_Dims_t
{
    bool  exists;
    float thick;
    float zdown;
    Aperture_Dims_t();
};

// RangeShifter_Dims_t ---------------------------------------
struct RangeShifter_Dims_t
{
    bool  exists;
    float thick;
    float zdown;
    float zup;
    RangeShifter_Dims_t();
private:
    float density;
public:
    float wepl;
};

// SAD_t ---------------------------------------
struct SAD_t
{
    SAD_t();
    SAD_t(std::string machine);
    float a;
    float b;
};


#endif
