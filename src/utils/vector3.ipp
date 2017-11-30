#include <cmath>
#include <iostream>
#include <vector>

#include "vector2.hpp"
#include "vector4.hpp"

// Vector3_t ---------------------------------------
template<class T>
Vector3_t<T>::Vector3_t()
{
    x = 0;
    y = 0;
    z = 0;
}

template<class T>
Vector3_t<T>::Vector3_t(T a)
{
    x = a;
    y = a;
    z = a;
}

template<class T>
Vector3_t<T>::Vector3_t(T a, T b, T c)
{
    x = a;
    y = b;
    z = c;
}

template<class T>
Vector3_t<T>::Vector3_t(const std::vector<T> &v)
{
    x = v.at(0);
    y = v.at(1);
    z = v.at(2);
}

template<class T>
Vector3_t<T>::Vector3_t(const Vector3_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
    z = obj.z;
}

template<class T>
Vector3_t<T>::Vector3_t(const Vector2_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
    z = 0;
}

template<class T>
Vector3_t<T>::Vector3_t(const Vector4_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
    z = obj.z;
}

template<class T>
void Vector3_t<T>::print(std::string sep) const
{
    std::cout << x << "  ";
    std::cout << y << "  ";
    std::cout << z << sep;
}

template<class T>
float Vector3_t<T>::length() const
{
    return std::sqrt(x*x + y*y + z*z);
}

template<class T>
float Vector3_t<T>::length2() const
{
    return x*x + y*y + z*z;
}

template<class T>
float Vector3_t<T>::dot(const Vector3_t<T>& a) const
{
    return x*a.x + y*a.y + z*a.z;
}

template<class T>
Vector3_t<T> Vector3_t<T>::cross (const Vector3_t<T>& v) const
{
    return Vector3_t<T>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x);
}

template<class T>
void Vector3_t<T>::normalize()
{
    float d = this->length();
    x /= d;
    y /= d;
    z /= d;
}

template<class T>
Vector3_t<T> Vector3_t<T>::get_normalized() const
{
    float d = this->length();
    return Vector3_t<T>(x/d, y/d, z/d);
}

template<class T>
void Vector3_t<T>::rotate(const float& gantry, const float& couch)
{
    float c_couch  = cos(couch);
    float s_couch  = sin(couch);
    float c_gantry = cos(gantry);
    float s_gantry = sin(gantry);

    Vector3_t<T> temp(x, y, z);
    x = temp.x*c_couch - s_couch*(temp.y*s_gantry + temp.z*c_gantry);
    y = temp.y*c_gantry - temp.z*s_gantry;
    z = temp.x*s_couch + c_couch*(temp.y*s_gantry + temp.z*c_gantry);
}

template<class T>
Vector3_t<T> Vector3_t<T>::get_rotated(const float& gantry, const float& couch) const
{
    float c_couch  = cos(couch);
    float s_couch  = sin(couch);
    float c_gantry = cos(gantry);
    float s_gantry = sin(gantry);

    Vector3_t<T> temp(x, y, z);
    float a_ = temp.x*c_couch - s_couch*(temp.y*s_gantry + temp.z*c_gantry);
    float b_ = temp.y*c_gantry - temp.z*s_gantry;
    float c_ = temp.x*s_couch + c_couch*(temp.y*s_gantry + temp.z*c_gantry);
    return Vector3_t<T>(a_, b_, c_);
}

// OPERATORS AND ACCESS
template<class T>
const T& Vector3_t<T>::at(int idx) const
{
    idx = idx % 3;
    if (idx == 0)
        return x;
    if (idx == 1)
        return y;
    if (idx == 2)
        return z;
    return x; // it will never get to this line
}

template<class T>
const T& Vector3_t<T>::operator [](int idx) const
{
    if (idx == 0)
        return x;
    if (idx == 1)
        return y;
    if (idx == 2)
        return z;
    std::cerr << "ERROR! Index out of range!" << std::endl;
    exit(EXIT_FAILURE);
    return x;
}

template<class T> template<class U>
Vector3_t<T>& Vector3_t<T>::operator+=(const Vector3_t<U>& rhs)
{
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
}
template<class T> template<typename U>
Vector3_t<T>& Vector3_t<T>::operator+=(const U& rhs)
{
    x += rhs;
    y += rhs;
    z += rhs;
    return *this;
}
template<class T> template<class U>
Vector3_t<T>& Vector3_t<T>::operator-=(const Vector3_t<U>& rhs)
{
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
}
template<class T> template<typename U>
Vector3_t<T>& Vector3_t<T>::operator-=(const U& rhs)
{
    x -= rhs;
    y -= rhs;
    z -= rhs;
    return *this;
}
template<class T> template<class U>
Vector3_t<T>& Vector3_t<T>::operator*=(const Vector3_t<U>& rhs)
{
    x *= rhs.x;
    y *= rhs.y;
    z *= rhs.z;
    return *this;
}
template<class T> template<typename U>
Vector3_t<T>& Vector3_t<T>::operator*=(const U& rhs)
{
    x *= rhs;
    y *= rhs;
    z *= rhs;
    return *this;
}
template<class T> template<class U>
Vector3_t<T>& Vector3_t<T>::operator/=(const Vector3_t<U>& rhs)
{
    x /= rhs.x;
    y /= rhs.y;
    z /= rhs.z;
    return *this;
}
template<class T> template<typename U>
Vector3_t<T>& Vector3_t<T>::operator/=(const U& rhs)
{
    x /= rhs;
    y /= rhs;
    z /= rhs;
    return *this;
}

template<class T, typename U>
Vector3_t<T> operator+(Vector3_t<T> lhs, const U& rhs)
{
    return lhs += rhs;
}
template<class T, typename U>
Vector3_t<T> operator+(const U& rhs, Vector3_t<T> lhs)
{
    return lhs += rhs;
}
template<class T, typename U>
Vector3_t<T> operator+(Vector3_t<T> lhs, const Vector3_t<U>& rhs)
{
    return lhs += rhs;
}

template<class T, typename U>
Vector3_t<T> operator-(Vector3_t<T> lhs, const U& rhs)
{
    return lhs -= rhs;
}
template<class T, typename U>
Vector3_t<T> operator-(const U& rhs, Vector3_t<T> lhs)
{
    return lhs -= rhs;
}
template<class T, typename U>
Vector3_t<T> operator-(Vector3_t<T> lhs, const Vector3_t<U>& rhs)
{
    return lhs -= rhs;
}

template<class T, typename U>
Vector3_t<T> operator*(Vector3_t<T> lhs, const U& rhs)
{
    return lhs *= rhs;
}
template<class T, typename U>
Vector3_t<T> operator*(const U& rhs, Vector3_t<T> lhs)
{
    return lhs *= rhs;
}
template<class T, typename U>
Vector3_t<T> operator*(Vector3_t<T> lhs, const Vector3_t<U>& rhs)
{
    return lhs *= rhs;
}

template<class T, typename U>
Vector3_t<T> operator/(Vector3_t<T> lhs, const U& rhs)
{
    return lhs /= rhs;
}
template<class T, typename U>
Vector3_t<T> operator/(const U& rhs, Vector3_t<T> lhs)
{
    return lhs /= rhs;
}
template<class T, typename U>
Vector3_t<T> operator/(Vector3_t<T> lhs, const Vector3_t<U>& rhs)
{
    return lhs /= rhs;
}
