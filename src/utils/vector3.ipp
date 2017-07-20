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
void Vector3_t<T>::print()
{
    std::cout << x << "\t";
    std::cout << y << "\t";
    std::cout << z << std::endl;
}

template<class T>
double Vector3_t<T>::length()
{
    return std::sqrt(x*x + y*y + z*z);
}

template<class T>
double Vector3_t<T>::length2()
{
    return x*x + y*y + z*z;
}

template<class T>
double Vector3_t<T>::dot(const Vector3_t<T>& a)
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
    double d = this->length();
    x /= d;
    y /= d;
    z /= d;
}

template<class T>
void Vector3_t<T>::rotate(const double& gantry, const double& couch)
{
    double c_couch  = cos(couch);
    double s_couch  = sin(couch);
    double c_gantry = cos(gantry);
    double s_gantry = sin(gantry);

    Vector3_t<T> temp(x, y, z);
    x = temp.x*c_couch - s_couch*(temp.y*s_gantry + temp.z*c_gantry);
    y = temp.y*c_gantry - temp.z*s_gantry;
    z = temp.x*s_couch + c_couch*(temp.y*s_gantry + temp.z*c_gantry);
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
