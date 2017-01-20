#ifndef __SPECIAL_TYPES_HPP__
#define __SPECIAL_TYPES_HPP__

#include <vector>

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

#endif