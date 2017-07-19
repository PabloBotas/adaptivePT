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
    Vector4_t(double a);
    Vector4_t(double a, double b, double c, double d);
    Vector4_t(const std::vector<T> &v);
    Vector4_t(const Vector4_t<T> &obj);
    Vector4_t(const Vector2_t<T> &obj);
    Vector4_t(const Vector3_t<T> &obj);
    Vector4_t operator+(const Vector4_t& rhs);
    Vector4_t operator-(const Vector4_t& rhs);
    Vector4_t operator*(const Vector4_t& rhs);
    Vector4_t operator/(const Vector4_t& rhs);
    void rotate(const double& gantry, const double& couch);
    const T& operator [](int idx) const;
    void print();

    T x;
    T y;
    T z;
    T w;
};

template <class T>
using Array4 = std::vector< Vector4_t<T> >;

inline Vector4_t<double> operator+(const Vector4_t<double>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}
inline Vector4_t<double> operator+(const Vector4_t<double>& a, const Vector4_t<int>& b)
{
    return Vector4_t<double>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}
inline Vector4_t<double> operator+(const Vector4_t<double>& a, const Vector4_t<unsigned int>& b)
{
    return Vector4_t<double>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}
inline Vector4_t<double> operator+(const Vector4_t<int>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}
inline Vector4_t<double> operator+(const Vector4_t<unsigned int>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}
inline Vector4_t<double> operator+(const Vector4_t<double>& a, const double& b)
{
    return Vector4_t<double>(a.x+b, a.y+b, a.z+b, a.w+b);
}
inline Vector4_t<double> operator+(const double& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a+b.x, a+b.y, a+b.z, a+b.w);
}

inline Vector4_t<double> operator-(const Vector4_t<double>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}
inline Vector4_t<double> operator-(const Vector4_t<double>& a, const Vector4_t<int>& b)
{
    return Vector4_t<double>(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}
inline Vector4_t<double> operator-(const Vector4_t<double>& a, const Vector4_t<unsigned int>& b)
{
    return Vector4_t<double>(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}
inline Vector4_t<double> operator-(const Vector4_t<int>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}
inline Vector4_t<double> operator-(const Vector4_t<unsigned int>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}
inline Vector4_t<double> operator-(const Vector4_t<double>& a, const double& b)
{
    return Vector4_t<double>(a.x-b, a.y-b, a.z-b, a.w-b);
}
inline Vector4_t<double> operator-(const double& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a-b.x, a-b.y, a-b.z, a-b.w);
}

inline Vector4_t<double> operator*(const Vector4_t<double>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}
inline Vector4_t<double> operator*(const Vector4_t<double>& a, const Vector4_t<int>& b)
{
    return Vector4_t<double>(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}
inline Vector4_t<double> operator*(const Vector4_t<double>& a, const Vector4_t<unsigned int>& b)
{
    return Vector4_t<double>(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}
inline Vector4_t<double> operator*(const Vector4_t<int>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}
inline Vector4_t<double> operator*(const Vector4_t<unsigned int>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}
inline Vector4_t<double> operator*(const Vector4_t<double>& a, const double& b)
{
    return Vector4_t<double>(a.x*b, a.y*b, a.z*b, a.w*b);
}
inline Vector4_t<double> operator*(const double& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a*b.x, a*b.y, a*b.z, a*b.w);
}

inline Vector4_t<double> operator/(const Vector4_t<double>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}
inline Vector4_t<double> operator/(const Vector4_t<double>& a, const Vector4_t<int>& b)
{
    return Vector4_t<double>(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}
inline Vector4_t<double> operator/(const Vector4_t<double>& a, const Vector4_t<unsigned int>& b)
{
    return Vector4_t<double>(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}
inline Vector4_t<double> operator/(const Vector4_t<int>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}
inline Vector4_t<double> operator/(const Vector4_t<unsigned int>& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}
inline Vector4_t<double> operator/(const Vector4_t<double>& a, const double& b)
{
    return Vector4_t<double>(a.x/b, a.y/b, a.z/b, a.w/b);
}
inline Vector4_t<double> operator/(const double& a, const Vector4_t<double>& b)
{
    return Vector4_t<double>(a/b.x, a/b.y, a/b.z, a/b.w);
}

inline Vector4_t<float> operator+(const Vector4_t<float>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}
inline Vector4_t<float> operator+(const Vector4_t<float>& a, const Vector4_t<int>& b)
{
    return Vector4_t<float>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}
inline Vector4_t<float> operator+(const Vector4_t<float>& a, const Vector4_t<unsigned int>& b)
{
    return Vector4_t<float>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}
inline Vector4_t<float> operator+(const Vector4_t<int>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}
inline Vector4_t<float> operator+(const Vector4_t<unsigned int>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}
inline Vector4_t<float> operator+(const Vector4_t<float>& a, const float& b)
{
    return Vector4_t<float>(a.x+b, a.y+b, a.z+b, a.w+b);
}
inline Vector4_t<float> operator+(const float& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a+b.x, a+b.y, a+b.z, a+b.w);
}

inline Vector4_t<float> operator-(const Vector4_t<float>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}
inline Vector4_t<float> operator-(const Vector4_t<float>& a, const Vector4_t<int>& b)
{
    return Vector4_t<float>(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}
inline Vector4_t<float> operator-(const Vector4_t<float>& a, const Vector4_t<unsigned int>& b)
{
    return Vector4_t<float>(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}
inline Vector4_t<float> operator-(const Vector4_t<int>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}
inline Vector4_t<float> operator-(const Vector4_t<unsigned int>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}
inline Vector4_t<float> operator-(const Vector4_t<float>& a, const float& b)
{
    return Vector4_t<float>(a.x-b, a.y-b, a.z-b, a.w-b);
}
inline Vector4_t<float> operator-(const float& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a-b.x, a-b.y, a-b.z, a-b.w);
}

inline Vector4_t<float> operator*(const Vector4_t<float>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}
inline Vector4_t<float> operator*(const Vector4_t<float>& a, const Vector4_t<int>& b)
{
    return Vector4_t<float>(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}
inline Vector4_t<float> operator*(const Vector4_t<float>& a, const Vector4_t<unsigned int>& b)
{
    return Vector4_t<float>(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}
inline Vector4_t<float> operator*(const Vector4_t<int>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}
inline Vector4_t<float> operator*(const Vector4_t<unsigned int>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}
inline Vector4_t<float> operator*(const Vector4_t<float>& a, const float& b)
{
    return Vector4_t<float>(a.x*b, a.y*b, a.z*b, a.w*b);
}
inline Vector4_t<float> operator*(const float& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a*b.x, a*b.y, a*b.z, a*b.w);
}

inline Vector4_t<float> operator/(const Vector4_t<float>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}
inline Vector4_t<float> operator/(const Vector4_t<float>& a, const Vector4_t<int>& b)
{
    return Vector4_t<float>(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}
inline Vector4_t<float> operator/(const Vector4_t<float>& a, const Vector4_t<unsigned int>& b)
{
    return Vector4_t<float>(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}
inline Vector4_t<float> operator/(const Vector4_t<int>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}
inline Vector4_t<float> operator/(const Vector4_t<unsigned int>& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}
inline Vector4_t<float> operator/(const Vector4_t<float>& a, const float& b)
{
    return Vector4_t<float>(a.x/b, a.y/b, a.z/b, a.w/b);
}
inline Vector4_t<float> operator/(const float& a, const Vector4_t<float>& b)
{
    return Vector4_t<float>(a/b.x, a/b.y, a/b.z, a/b.w);
}

#endif
