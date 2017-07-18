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

// SAD_t ---------------------------------------
SAD_t::SAD_t()
{
    a = 0;
    b = 0;
}

SAD_t::SAD_t(std::string machine)
{
    if(machine.compare(utils::toLower("TopasSmallSpots"))      == 0 ||
       machine.compare(utils::toLower("Topas_a5_SmallSpots"))  == 0 ||
       machine.compare(utils::toLower("TopasMediumSpots"))     == 0 ||
       machine.compare(utils::toLower("Topas_a5_MediumSpots")) == 0)
    {
        a = 1940.f;
        b = 2340.f;
    }
    else if(machine.compare(utils::toLower("TopasMGHR4")) == 0)
    {
        a = 2340.f;
        b = 1940.f;
    }
}


