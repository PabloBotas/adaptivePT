#ifndef __SPECIAL_TYPES_HPP__
#define __SPECIAL_TYPES_HPP__

#include <vector>
#include <string>

template <class T>
class Vector3_t
{
public:
    Vector3_t();
    Vector3_t(const std::vector<T> &v);
    Vector3_t(const Vector3_t<T> &obj);

    T x;
    T y;
    T z;
};

template <class T>
class Vector4_t
{
public:
    Vector4_t();
    Vector4_t(const std::vector<T> &v);
    Vector4_t(const Vector4_t<T> &obj);
    Vector4_t(const Vector3_t<T> &obj);

    T x;
    T y;
    T z;
    T w;
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
    Vector3_t<float> offset;
    Vector3_t<float> isocenter;
    Vector3_t<float> d;
    Vector3_t<unsigned int> n;
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
