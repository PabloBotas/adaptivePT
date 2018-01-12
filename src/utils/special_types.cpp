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
                     adapted(false),
                     thick(0),
                     zdown(0),
                     zup(0),
                     density(1.15),
                     wepl(thick*density)
{
}

void RangeShifter_Dims_t::substract(float x)
{
    adapted = true;
    wepl -= x;
    density = wepl/thick;
    if (wepl == 0)
        exists = false;
}

void RangeShifter_Dims_t::create(float zup_, float w)
{
    // __RANGE_SHIFTER_THICKNESS__ is defined in CMake step, defaulting to 0.01
    exists = true;
    zup = zup_;
    zdown = zup - __RANGE_SHIFTER_THICKNESS__;
    thick = __RANGE_SHIFTER_THICKNESS__;
    wepl = w;
    density = wepl/thick;
}

void RangeShifter_Dims_t::add(float w)
{
    adapted = true;
    wepl += w;
    density = w/thick;
}

void RangeShifter_Dims_t::set_adapted()
{
    adapted = true;
}

void RangeShifter_Dims_t::set_wepl(float w)
{
    adapted = true;
    wepl = w;
    density = w/thick;
}

std::string RangeShifter_Dims_t::get_info_as_str()
{
    std::string str;
    if (exists) {
        str = "adapted " + std::string(adapted ? "true" : "false") + " ";
        str += "thick " + std::to_string(thick) + " ";
        str += "wepl " + std::to_string(wepl) + " ";
        str += "density " + std::to_string(density) + " ";
        str += "zdown " + std::to_string(zdown) + " ";
        str += "zup " + std::to_string(zup);
    } else {
        str = "Empty";
    }
    return str;
}

void RangeShifter_Dims_t::print(std::ostream& stream)
{
    stream << get_info_as_str() << std::endl;
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


