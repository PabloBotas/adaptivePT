#ifndef __VECTOR2_HPP__
#define __VECTOR2_HPP__

#include "vector3.hpp"
#include "vector4.hpp"
#include <vector>

template <class T> class Vector3_t;
template <class T> class Vector4_t;

template <class T>
class Vector2_t
{
public:
    Vector2_t();
    Vector2_t(double a, double b);
    Vector2_t(const std::vector<T> &v);
    Vector2_t(const Vector2_t<T> &obj);
    Vector2_t(const Vector3_t<T> &obj);
    Vector2_t(const Vector4_t<T> &obj);
    void print();

    T x;
    T y;
};

template <class T>
using Array2 = std::vector< Vector2_t<T> >;

inline Vector2_t<double> operator+(const Vector2_t<double>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x+b.x, a.y+b.y);
}
inline Vector2_t<double> operator+(const Vector2_t<double>& a, const Vector2_t<int>& b)
{
    return Vector2_t<double>(a.x+b.x, a.y+b.y);
}
inline Vector2_t<double> operator+(const Vector2_t<double>& a, const Vector2_t<unsigned int>& b)
{
    return Vector2_t<double>(a.x+b.x, a.y+b.y);
}
inline Vector2_t<double> operator+(const Vector2_t<int>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x+b.x, a.y+b.y);
}
inline Vector2_t<double> operator+(const Vector2_t<unsigned int>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x+b.x, a.y+b.y);
}
inline Vector2_t<double> operator+(const Vector2_t<double>& a, const double& b)
{
    return Vector2_t<double>(a.x+b, a.y+b);
}
inline Vector2_t<double> operator+(const double& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a+b.x, a+b.y);
}

inline Vector2_t<double> operator-(const Vector2_t<double>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x-b.x, a.y-b.y);
}
inline Vector2_t<double> operator-(const Vector2_t<double>& a, const Vector2_t<int>& b)
{
    return Vector2_t<double>(a.x-b.x, a.y-b.y);
}
inline Vector2_t<double> operator-(const Vector2_t<double>& a, const Vector2_t<unsigned int>& b)
{
    return Vector2_t<double>(a.x-b.x, a.y-b.y);
}
inline Vector2_t<double> operator-(const Vector2_t<int>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x-b.x, a.y-b.y);
}
inline Vector2_t<double> operator-(const Vector2_t<unsigned int>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x-b.x, a.y-b.y);
}
inline Vector2_t<double> operator-(const Vector2_t<double>& a, const double& b)
{
    return Vector2_t<double>(a.x-b, a.y-b);
}
inline Vector2_t<double> operator-(const double& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a-b.x, a-b.y);
}

inline Vector2_t<double> operator*(const Vector2_t<double>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x*b.x, a.y*b.y);
}
inline Vector2_t<double> operator*(const Vector2_t<double>& a, const Vector2_t<int>& b)
{
    return Vector2_t<double>(a.x*b.x, a.y*b.y);
}
inline Vector2_t<double> operator*(const Vector2_t<double>& a, const Vector2_t<unsigned int>& b)
{
    return Vector2_t<double>(a.x*b.x, a.y*b.y);
}
inline Vector2_t<double> operator*(const Vector2_t<int>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x*b.x, a.y*b.y);
}
inline Vector2_t<double> operator*(const Vector2_t<unsigned int>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x*b.x, a.y*b.y);
}
inline Vector2_t<double> operator*(const Vector2_t<double>& a, const double& b)
{
    return Vector2_t<double>(a.x*b, a.y*b);
}
inline Vector2_t<double> operator*(const double& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a*b.x, a*b.y);
}

inline Vector2_t<double> operator/(const Vector2_t<double>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x/b.x, a.y/b.y);
}
inline Vector2_t<double> operator/(const Vector2_t<double>& a, const Vector2_t<int>& b)
{
    return Vector2_t<double>(a.x/b.x, a.y/b.y);
}
inline Vector2_t<double> operator/(const Vector2_t<double>& a, const Vector2_t<unsigned int>& b)
{
    return Vector2_t<double>(a.x/b.x, a.y/b.y);
}
inline Vector2_t<double> operator/(const Vector2_t<int>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x/b.x, a.y/b.y);
}
inline Vector2_t<double> operator/(const Vector2_t<unsigned int>& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a.x/b.x, a.y/b.y);
}
inline Vector2_t<double> operator/(const Vector2_t<double>& a, const double& b)
{
    return Vector2_t<double>(a.x/b, a.y/b);
}
inline Vector2_t<double> operator/(const double& a, const Vector2_t<double>& b)
{
    return Vector2_t<double>(a/b.x, a/b.y);
}

#endif
