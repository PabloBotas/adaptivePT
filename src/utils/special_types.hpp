#ifndef __SPECIAL_TYPES_HPP__
#define __SPECIAL_TYPES_HPP__

#include "vector3.hpp"
#include <iostream>
#include <string>

// BeamAngles_t ---------------------------------------
struct BeamAngles_t
{
    float gantry;
    float couch;
};

// Planes_t ---------------------------------------
struct Planes_t
{
    Array3<float> dir;
    Array3<float> p;
    Array3<float> source_a;
    Array3<float> source_b;
    Planes_t();
    Planes_t(size_t n);
};

// CT_Dims_t ---------------------------------------
struct CT_Dims_t
{
    Vector3_t<float> offset;
    Vector3_t<float> origin;
    Vector3_t<float> isocenter;
    Vector3_t<float> d;
    Vector3_t<unsigned int> n;
    unsigned int total;
};

// Aperture_Dims_t ---------------------------------------
struct Aperture_Dims_t
{
    bool exists;
    float thick;
    float zdown;
    Aperture_Dims_t();
};

// RangeShifter_Dims_t ---------------------------------------
struct RangeShifter_Dims_t
{
    bool exists;
    bool adapted;
    float thick;
    float zdown;
    float zup;
    float density;
    RangeShifter_Dims_t();
    float wepl;
    void add (float w);
    void set_adapted();
    void set_wepl(float w);
    void substract (float w);
    void create (float pos, float w);
    std::string get_info_as_str();
    void print(std::ostream& stream = std::cout);
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
