#include "vector4.hpp"

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
Vector4_t<T>::Vector4_t(const std::vector<T> &v)
{
    x = v.at(0);
    y = v.at(1);
    z = v.at(2);
    z = v.at(4);
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

template class Vector4_t<float>;
