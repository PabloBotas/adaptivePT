#ifndef __VECTOR4_HPP__
#define __VECTOR4_HPP__

// #include "vector2.hpp"
// #include "vector3.hpp"
#include <vector>
#include <string>

template <class T> class Vector2_t;
template <class T> class Vector3_t;

template <class T>
class Vector4_t
{
public:
    Vector4_t();
    Vector4_t(T a);
    Vector4_t(T a, T b, T c, T d);
    Vector4_t(const std::vector<T> &v);
    Vector4_t(const Vector4_t<T> &obj);
    Vector4_t(const Vector2_t<T> &obj);
    Vector4_t(const Vector3_t<T> &obj);
    void print(std::string sep = "\n") const;
    void print_as_3D(std::string sep = "\n") const;
    double length() const;
    double length2() const;
    double dot(const Vector4_t<T>& a) const;
    Vector4_t<T> cross(const Vector4_t<T>& v) const;
    void normalize();
    Vector4_t<T> get_normalized() const;
    void rotate(const double& gantry, const double& couch);
    Vector4_t<T> get_rotated(const double& gantry, const double& couch) const;

    T x;
    T y;
    T z;
    T w;

    const T& operator [](int idx) const;
    const T& at(int idx) const;
    template <class U> Vector4_t<T>& operator+=(const Vector4_t<U>& rhs);
    template <class U> Vector4_t<T>& operator+=(const Vector3_t<U>& rhs);
    template <typename U> Vector4_t<T>& operator+=(const U& rhs);
    template <class U> Vector4_t<T>& operator-=(const Vector4_t<U>& rhs);
    template <class U> Vector4_t<T>& operator-=(const Vector3_t<U>& rhs);
    template <typename U> Vector4_t<T>& operator-=(const U& rhs);
    template <class U> Vector4_t<T>& operator*=(const Vector4_t<U>& rhs);
    template <class U> Vector4_t<T>& operator*=(const Vector3_t<U>& rhs);
    template <typename U> Vector4_t<T>& operator*=(const U& rhs);
    template <class U> Vector4_t<T>& operator/=(const Vector4_t<U>& rhs);
    template <class U> Vector4_t<T>& operator/=(const Vector3_t<U>& rhs);
    template <typename U> Vector4_t<T>& operator/=(const U& rhs);
};

template <class T>
using Array4 = std::vector< Vector4_t<T> >;

#include "vector4.ipp"

#endif
