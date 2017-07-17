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
    float length();
    void normalize();

    T x;
    T y;
    T z;
};

template <class T>
using Array3 = std::vector< Vector3_t<T> >;

inline Vector3_t<float> operator+(const Vector3_t<float>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x+b.x, a.y+b.y, a.z+b.z);
}
inline Vector3_t<float> operator+(const Vector3_t<float>& a, const Vector3_t<int>& b)
{
    return Vector3_t<float>(a.x+b.x, a.y+b.y, a.z+b.z);
}
inline Vector3_t<float> operator+(const Vector3_t<float>& a, const Vector3_t<unsigned int>& b)
{
    return Vector3_t<float>(a.x+b.x, a.y+b.y, a.z+b.z);
}
inline Vector3_t<float> operator+(const Vector3_t<int>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x+b.x, a.y+b.y, a.z+b.z);
}
inline Vector3_t<float> operator+(const Vector3_t<unsigned int>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x+b.x, a.y+b.y, a.z+b.z);
}
inline Vector3_t<float> operator+(const Vector3_t<float>& a, const float& b)
{
    return Vector3_t<float>(a.x+b, a.y+b, a.z+b);
}
inline Vector3_t<float> operator+(const float& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a+b.x, a+b.y, a+b.z);
}

inline Vector3_t<float> operator-(const Vector3_t<float>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x-b.x, a.y-b.y, a.z-b.z);
}
inline Vector3_t<float> operator-(const Vector3_t<float>& a, const Vector3_t<int>& b)
{
    return Vector3_t<float>(a.x-b.x, a.y-b.y, a.z-b.z);
}
inline Vector3_t<float> operator-(const Vector3_t<float>& a, const Vector3_t<unsigned int>& b)
{
    return Vector3_t<float>(a.x-b.x, a.y-b.y, a.z-b.z);
}
inline Vector3_t<float> operator-(const Vector3_t<int>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x-b.x, a.y-b.y, a.z-b.z);
}
inline Vector3_t<float> operator-(const Vector3_t<unsigned int>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x-b.x, a.y-b.y, a.z-b.z);
}
inline Vector3_t<float> operator-(const Vector3_t<float>& a, const float& b)
{
    return Vector3_t<float>(a.x-b, a.y-b, a.z-b);
}
inline Vector3_t<float> operator-(const float& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a-b.x, a-b.y, a-b.z);
}

inline Vector3_t<float> operator*(const Vector3_t<float>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x*b.x, a.y*b.y, a.z*b.z);
}
inline Vector3_t<float> operator*(const Vector3_t<float>& a, const Vector3_t<int>& b)
{
    return Vector3_t<float>(a.x*b.x, a.y*b.y, a.z*b.z);
}
inline Vector3_t<float> operator*(const Vector3_t<float>& a, const Vector3_t<unsigned int>& b)
{
    return Vector3_t<float>(a.x*b.x, a.y*b.y, a.z*b.z);
}
inline Vector3_t<float> operator*(const Vector3_t<int>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x*b.x, a.y*b.y, a.z*b.z);
}
inline Vector3_t<float> operator*(const Vector3_t<unsigned int>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x*b.x, a.y*b.y, a.z*b.z);
}
inline Vector3_t<float> operator*(const Vector3_t<float>& a, const float& b)
{
    return Vector3_t<float>(a.x*b, a.y*b, a.z*b);
}
inline Vector3_t<float> operator*(const float& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a*b.x, a*b.y, a*b.z);
}

inline Vector3_t<float> operator/(const Vector3_t<float>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x/b.x, a.y/b.y, a.z/b.z);
}
inline Vector3_t<float> operator/(const Vector3_t<float>& a, const Vector3_t<int>& b)
{
    return Vector3_t<float>(a.x/b.x, a.y/b.y, a.z/b.z);
}
inline Vector3_t<float> operator/(const Vector3_t<float>& a, const Vector3_t<unsigned int>& b)
{
    return Vector3_t<float>(a.x/b.x, a.y/b.y, a.z/b.z);
}
inline Vector3_t<float> operator/(const Vector3_t<int>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x/b.x, a.y/b.y, a.z/b.z);
}
inline Vector3_t<float> operator/(const Vector3_t<unsigned int>& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a.x/b.x, a.y/b.y, a.z/b.z);
}
inline Vector3_t<float> operator/(const Vector3_t<float>& a, const float& b)
{
    return Vector3_t<float>(a.x/b, a.y/b, a.z/b);
}
inline Vector3_t<float> operator/(const float& a, const Vector3_t<float>& b)
{
    return Vector3_t<float>(a/b.x, a/b.y, a/b.z);
}

#endif
