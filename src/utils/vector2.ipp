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
Vector2_t<T>::Vector2_t(T a)
{
    x = a;
    y = a;
}

template<class T>
Vector2_t<T>::Vector2_t(T a, T b)
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
void Vector2_t<T>::print(std::string sep) const
{
        std::cout << x << "  ";
        std::cout << y << sep;
}


// OPERATORS AND ACCESS
template<class T>
const T& Vector2_t<T>::at(int idx) const
{
    idx = idx % 2;
    if (idx == 0)
        return x;
    if (idx == 1)
        return y;
    return x; // it will never get to this line
}

template<class T>
const T& Vector2_t<T>::operator [](int idx) const
{
    if (idx == 0)
        return x;
    if (idx == 1)
        return y;
    std::cerr << "ERROR! Index out of range!" << std::endl;
    exit(EXIT_FAILURE);
    return x;
}

template<class T> template<class U>
Vector2_t<T>& Vector2_t<T>::operator+=(const Vector2_t<U>& rhs)
{
    x += rhs.x;
    y += rhs.y;
    return *this;
}
template<class T> template<typename U>
Vector2_t<T>& Vector2_t<T>::operator+=(const U& rhs)
{
    x += rhs;
    y += rhs;
    return *this;
}
template<class T> template<class U>
Vector2_t<T>& Vector2_t<T>::operator-=(const Vector2_t<U>& rhs)
{
    x -= rhs.x;
    y -= rhs.y;
    return *this;
}
template<class T> template<typename U>
Vector2_t<T>& Vector2_t<T>::operator-=(const U& rhs)
{
    x -= rhs;
    y -= rhs;
    return *this;
}
template<class T> template<class U>
Vector2_t<T>& Vector2_t<T>::operator*=(const Vector2_t<U>& rhs)
{
    x *= rhs.x;
    y *= rhs.y;
    return *this;
}
template<class T> template<typename U>
Vector2_t<T>& Vector2_t<T>::operator*=(const U& rhs)
{
    x *= rhs;
    y *= rhs;
    return *this;
}
template<class T> template<class U>
Vector2_t<T>& Vector2_t<T>::operator/=(const Vector2_t<U>& rhs)
{
    x /= rhs.x;
    y /= rhs.y;
    return *this;
}
template<class T> template<typename U>
Vector2_t<T>& Vector2_t<T>::operator/=(const U& rhs)
{
    x /= rhs;
    y /= rhs;
    return *this;
}

template<class T, typename U>
Vector2_t<T> operator+(Vector2_t<T> lhs, const U& rhs)
{
    return lhs += rhs;
}
template<class T, typename U>
Vector2_t<T> operator+(const U& rhs, Vector2_t<T> lhs)
{
    return lhs += rhs;
}
template<class T, typename U>
Vector2_t<T> operator+(Vector2_t<T> lhs, const Vector2_t<U>& rhs)
{
    return lhs += rhs;
}

template<class T, typename U>
Vector2_t<T> operator-(Vector2_t<T> lhs, const U& rhs)
{
    return lhs -= rhs;
}
template<class T, typename U>
Vector2_t<T> operator-(const U& rhs, Vector2_t<T> lhs)
{
    return lhs -= rhs;
}
template<class T, typename U>
Vector2_t<T> operator-(Vector2_t<T> lhs, const Vector2_t<U>& rhs)
{
    return lhs -= rhs;
}

template<class T, typename U>
Vector2_t<T> operator*(Vector2_t<T> lhs, const U& rhs)
{
    return lhs *= rhs;
}
template<class T, typename U>
Vector2_t<T> operator*(const U& rhs, Vector2_t<T> lhs)
{
    return lhs *= rhs;
}
template<class T, typename U>
Vector2_t<T> operator*(Vector2_t<T> lhs, const Vector2_t<U>& rhs)
{
    return lhs *= rhs;
}

template<class T, typename U>
Vector2_t<T> operator/(Vector2_t<T> lhs, const U& rhs)
{
    return lhs /= rhs;
}
template<class T, typename U>
Vector2_t<T> operator/(const U& rhs, Vector2_t<T> lhs)
{
    return lhs /= rhs;
}
template<class T, typename U>
Vector2_t<T> operator/(Vector2_t<T> lhs, const Vector2_t<U>& rhs)
{
    return lhs /= rhs;
}

