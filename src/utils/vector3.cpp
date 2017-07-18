#include "vector3.hpp"

#include <cmath>
#include <iostream>
#include <vector>

// Vector3_t ---------------------------------------
template<class T>
Vector3_t<T>::Vector3_t()
{
    x = 0;
    y = 0;
    z = 0;
}

template<class T>
Vector3_t<T>::Vector3_t(float a, float b, float c)
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
void Vector3_t<T>::print()
{
    std::cout << x << "\t";
    std::cout << y << "\t";
    std::cout << z << std::endl;
}

template<class T>
float Vector3_t<T>::length()
{
    return std::sqrt(x*x + y*y + z*z);
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
const T& Vector3_t<T>::operator [](int idx) const
{
    if (idx < 0)
        idx = 3 + (idx % 3);
    
    idx = idx % 3;
    if (idx == 0)
        return x;
    if (idx == 1)
        return y;
    if (idx == 2)
        return z;
    return x; // it will never get to this line
}

template class Vector3_t<int>;
template class Vector3_t<unsigned int>;
template class Vector3_t<float>;

