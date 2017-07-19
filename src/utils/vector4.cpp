#include "vector4.hpp"

#include <iostream>
#include <vector>
#include <cmath>

// Vector4_t ---------------------------------------
template<class T>
Vector4_t<T>::Vector4_t()
{
    x = 0;
    y = 0;
    z = 0;
    w = 0;
}

template<class T>
Vector4_t<T>::Vector4_t(double a)
{
    x = a;
    y = a;
    z = a;
    w = a;
}

template<class T>
Vector4_t<T>::Vector4_t(double a, double b, double c, double d)
{
    x = a;
    y = b;
    z = c;
    w = d;
}

template<class T>
Vector4_t<T>::Vector4_t(const std::vector<T> &v)
{
    x = v.at(0);
    y = v.at(1);
    z = v.at(2);
    w = v.at(4);
}

template<class T>
Vector4_t<T>::Vector4_t(const Vector4_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
    z = obj.z;
    w = obj.w;
}

template<class T>
Vector4_t<T>::Vector4_t(const Vector2_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
    z = 0;
    w = 0;
}

template<class T>
Vector4_t<T>::Vector4_t(const Vector3_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
    z = obj.z;
    w = 0;
}

template<class T>
void Vector4_t<T>::print()
{
    std::cout << x << "\t";
    std::cout << y << "\t";
    std::cout << z << "\t";
    std::cout << w << std::endl;
}

template<class T>
void Vector4_t<T>::rotate(const double& gantry, const double& couch)
{
    double c_couch  = cos(couch);
    double s_couch  = sin(couch);
    double c_gantry = cos(gantry);
    double s_gantry = sin(gantry);

    Vector4_t<T> temp(x, y, z, w);
    x = temp.x*c_couch - s_couch*(temp.y*s_gantry + temp.z*c_gantry);
    y = temp.y*c_gantry - temp.z*s_gantry;
    z = temp.x*s_couch + c_couch*(temp.y*s_gantry + temp.z*c_gantry);
    w = temp.w;
}

template<class T>
const T& Vector4_t<T>::operator [](int idx) const
{
    if (idx < 0)
        idx = 4 + (idx % 4);
    
    idx = idx % 4;
    if (idx == 0)
        return x;
    if (idx == 1)
        return y;
    if (idx == 2)
        return z;
    if (idx == 3)
        return w;
    return x; // it will never get to this line
}

template class Vector4_t<float>;
template class Vector4_t<double>;
