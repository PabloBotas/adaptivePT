#include "special_types.hpp"
#include "utils.hpp"

#include <vector>
#include <algorithm>

// Planes_t ---------------------------------------
Planes_t::Planes_t(){}

Planes_t::Planes_t(size_t n)
{
    dir.resize(n);
    p.resize(n);
    source_a.resize(n);
    source_b.resize(n);
}

// Aperture_Dims_t ---------------------------------------
Aperture_Dims_t::Aperture_Dims_t():
                 exists(false),
                 thick(0),
                 zdown(0)
{
}

// RangeShifter_Dims_t ---------------------------------------
RangeShifter_Dims_t::RangeShifter_Dims_t():
                     exists(false),
                     thick(0),
                     zdown(0),
                     zup(0),
                     density(1.15),
                     wepl(thick*density)
{
}

void RangeShifter_Dims_t::substract(float x)
{
    wepl -= x;
    density = wepl/thick;
    if (wepl == 0)
        exists = false;
}

void RangeShifter_Dims_t::create(float pos, float w)
{
    // __RANGE_SHIFTER_THICKNESS__ is defined in CMake step, defaulting to 0.01
    exists = true;
    zup = pos;
    zdown = pos - __RANGE_SHIFTER_THICKNESS__;
    thick = zup - zdown;
    density = w/thick;
    wepl = w;
}

void RangeShifter_Dims_t::add(float w)
{
    wepl += w;
    density = w/thick;
}

void RangeShifter_Dims_t::set_wepl(float w)
{
    wepl = w;
    density = w/thick;
}

// SAD_t ---------------------------------------
SAD_t::SAD_t()
{
    a = 0;
    b = 0;
}

SAD_t::SAD_t(std::string machine)
{
    if (machine == utils::toLower("TopasSmallSpots")     ||
        machine == utils::toLower("Topas_a5_SmallSpots") ||
        machine == utils::toLower("TopasMediumSpots")    ||
        machine == utils::toLower("Topas_a5_MediumSpots")) {
        a = 1940.f;
        b = 2340.f;
    } else if (machine == utils::toLower("TopasMGHR4") ||
               machine == utils::toLower("TopasMGHR5")) {
        a = 2340.f;
        b = 1940.f;
    }
}


