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
    Vector2_t(T a);
    Vector2_t(T a, T b);
    Vector2_t(const std::vector<T> &v);
    Vector2_t(const Vector2_t<T> &obj);
    Vector2_t(const Vector3_t<T> &obj);
    Vector2_t(const Vector4_t<T> &obj);
    void print(std::string sep = "\n") const;

    T x;
    T y;

    const T& operator [](int idx) const;
    const T& at(int idx) const;
    template <class U> Vector2_t<T>& operator+=(const Vector2_t<U>& rhs);
    template <typename U> Vector2_t<T>& operator+=(const U& rhs);
    template <class U> Vector2_t<T>& operator-=(const Vector2_t<U>& rhs);
    template <typename U> Vector2_t<T>& operator-=(const U& rhs);
    template <class U> Vector2_t<T>& operator*=(const Vector2_t<U>& rhs);
    template <typename U> Vector2_t<T>& operator*=(const U& rhs);
    template <class U> Vector2_t<T>& operator/=(const Vector2_t<U>& rhs);
    template <typename U> Vector2_t<T>& operator/=(const U& rhs);
};

template <class T>
using Array2 = std::vector< Vector2_t<T> >;

#include "vector2.ipp"

#endif
