#ifndef __SPECIAL_TYPES_HPP__
#define __SPECIAL_TYPES_HPP__

#include "vector3.hpp"
#include <string>

// BeamAngles_t ---------------------------------------
struct BeamAngles_t
{
    double gantry;
    double couch;
};

// Planes_t ---------------------------------------
struct Planes_t
{
    Array3<double> dir;
    Array3<double> p;
    Array3<double> source_a;
    Array3<double> source_b;
    Planes_t();
    Planes_t(size_t n);
};

// CT_Dims_t ---------------------------------------
struct CT_Dims_t
{
    Vector3_t<double> offset;
    Vector3_t<double> isocenter;
    Vector3_t<double> d;
    Vector3_t<unsigned int> n;
    unsigned int total;
};

// Aperture_Dims_t ---------------------------------------
struct Aperture_Dims_t
{
    bool exists;
    double thick;
    double zdown;
    Aperture_Dims_t();
};

// RangeShifter_Dims_t ---------------------------------------
struct RangeShifter_Dims_t
{
    bool exists;
    double thick;
    double zdown;
    double zup;
    RangeShifter_Dims_t();
private:
    double density;
public:
    double wepl;
};

// SAD_t ---------------------------------------
struct SAD_t
{
    SAD_t();
    SAD_t(std::string machine);
    double a;
    double b;
};


#endif
