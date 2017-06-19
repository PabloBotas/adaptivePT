#ifndef __VECTOR3_HPP__
#define __VECTOR3_HPP__

#include "vector2.hpp"
#include "vector4.hpp"
#include <vector>

template <class T> class Vector2_t;
template <class T> class Vector4_t;

template <class T>
class Vector3_t
{
public:
    Vector3_t();
    Vector3_t(float a, float b, float c);
    Vector3_t(const std::vector<T> &v);
    Vector3_t(const Vector3_t<T> &obj);
    Vector3_t(const Vector2_t<T> &obj);
    Vector3_t(const Vector4_t<T> &obj);
    void print();

    T x;
    T y;
    T z;
};

template <class T>
using Array3 = std::vector< Vector3_t<T> >;

#endif
