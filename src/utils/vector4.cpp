#include "vector4.hpp"

#include <iostream>
#include <vector>

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
Vector4_t<T>::Vector4_t(float a, float b, float c, float d)
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

template class Vector4_t<float>;
