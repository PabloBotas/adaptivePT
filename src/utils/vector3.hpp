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
    Vector3_t(double a);
    Vector3_t(double a, double b, double c);
    Vector3_t(const std::vector<T> &v);
    Vector3_t(const Vector3_t<T> &obj);
    Vector3_t(const Vector2_t<T> &obj);
    Vector3_t(const Vector4_t<T> &obj);
    void print();
    double length();
    double length2();
    double dot(const Vector3_t<T>& a);
    Vector3_t<T> cross(const Vector3_t<T>& v) const;
    void normalize();
    void rotate(const double& gantry, const double& couch);

    T x;
    T y;
    T z;

    const T& operator [](int idx) const;
    const T& at(int idx) const;
    template <class U> Vector3_t<T>& operator+=(const Vector3_t<U>& rhs);
    template <typename U> Vector3_t<T>& operator+=(const U& rhs);
    template <class U> Vector3_t<T>& operator-=(const Vector3_t<U>& rhs);
    template <typename U> Vector3_t<T>& operator-=(const U& rhs);
    template <class U> Vector3_t<T>& operator*=(const Vector3_t<U>& rhs);
    template <typename U> Vector3_t<T>& operator*=(const U& rhs);
    template <class U> Vector3_t<T>& operator/=(const Vector3_t<U>& rhs);
    template <typename U> Vector3_t<T>& operator/=(const U& rhs);
};

template <class T>
using Array3 = std::vector< Vector3_t<T> >;

#include "vector3.ipp"

#endif
