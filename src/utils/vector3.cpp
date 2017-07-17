#include "vector3.hpp"

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
    return srqt(x*x + y*y + z*z);
}

template<class T>
void Vector3_t<T>::normalize()
{
    float d = this->length();
    x = x/d;
    y = y/d;
    z = z/d;
}


template class Vector3_t<int>;
template class Vector3_t<unsigned int>;
template class Vector3_t<float>;
