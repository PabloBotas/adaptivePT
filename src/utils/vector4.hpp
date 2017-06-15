#ifndef __VECTOR4_HPP__
#define __VECTOR4_HPP__

#include "vector2.hpp"
#include "vector3.hpp"
#include <vector>

template <class T> class Vector2_t;
template <class T> class Vector3_t;

template <class T>
class Vector4_t
{
public:
    Vector4_t();
    Vector4_t(float a);
    Vector4_t(float a, float b, float c, float d);
    Vector4_t(const std::vector<T> &v);
    Vector4_t(const Vector4_t<T> &obj);
    Vector4_t(const Vector2_t<T> &obj);
    Vector4_t(const Vector3_t<T> &obj);
    void print();

    T x;
    T y;
    T z;
    T w;
};

#endif
