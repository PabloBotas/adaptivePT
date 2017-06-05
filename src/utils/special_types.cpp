#include "special_types.hpp"
#include "utils.hpp"

#include <vector>
#include <algorithm>

// Vector3_t ---------------------------------------
template<class T>
Vector3_t<T>::Vector3_t()
{
    x = 0;
    y = 0;
    z = 0;
}

template<class T>
Vector3_t<T>::Vector3_t(const std::vector<T> &v)
{
    x = v.at(0);
    y = v.at(1);
    z = v.at(2);
}

template<class T>
Vector3_t<T>::Vector3_t(const Vector3_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
    z = obj.z;
}

// Vector4_t ---------------------------------------
template<class T>
Vector4_t<T>::Vector4_t()
{
    x = 0;
    y = 0;
    z = 0;
    w = 0;
}

template<class T>
Vector4_t<T>::Vector4_t(const std::vector<T> &v)
{
    x = v.at(0);
    y = v.at(1);
    z = v.at(2);
    z = v.at(4);
}

template<class T>
Vector4_t<T>::Vector4_t(const Vector4_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
    z = obj.z;
    w = obj.w;
}

template<class T>
Vector4_t<T>::Vector4_t(const Vector3_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
    z = obj.z;
    w = 0;
}

template class Vector3_t<int>;
template class Vector3_t<unsigned int>;
template class Vector3_t<float>;
template class Vector4_t<float>;


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


