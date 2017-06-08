#include "vector2.hpp"

#include <iostream>
#include <vector>

// Vector2_t ---------------------------------------
template<class T>
Vector2_t<T>::Vector2_t()
{
    x = 0;
    y = 0;
}

template<class T>
Vector2_t<T>::Vector2_t(float a, float b)
{
    x = a;
    y = b;
}

template<class T>
Vector2_t<T>::Vector2_t(const std::vector<T> &v)
{
    x = v.at(0);
    y = v.at(1);
}

template<class T>
Vector2_t<T>::Vector2_t(const Vector2_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
}

template<class T>
Vector2_t<T>::Vector2_t(const Vector3_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
}

template<class T>
Vector2_t<T>::Vector2_t(const Vector4_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
}

template<class T>
void Vector2_t<T>::print()
{
        std::cout << x << "\t";
        std::cout << y << std::endl;
}


template class Vector2_t<short>;
