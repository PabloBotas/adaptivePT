#ifndef __SPECIAL_TYPES_HPP__
#define __SPECIAL_TYPES_HPP__

#include <vector>

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

#endif