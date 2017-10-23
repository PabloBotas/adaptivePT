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
Vector4_t<T>::Vector4_t(T a)
{
    x = a;
    y = a;
    z = a;
    w = a;
}

template<class T>
Vector4_t<T>::Vector4_t(T a, T b, T c, T d)
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
void Vector4_t<T>::print(std::string sep) const
{
    std::cout << x << "  ";
    std::cout << y << "  ";
    std::cout << z << "  ";
    std::cout << w << sep;
}

template<class T>
void Vector4_t<T>::print_as_3D(std::string sep) const
{
    std::cout << x << "  ";
    std::cout << y << "  ";
    std::cout << z << sep;
}

template<class T>
double Vector4_t<T>::length() const
{
    return std::sqrt(x*x + y*y + z*z);
}

template<class T>
double Vector4_t<T>::length2() const
{
    return x*x + y*y + z*z;
}

template<class T>
double Vector4_t<T>::dot(const Vector4_t<T>& a) const
{
    return x*a.x + y*a.y + z*a.z;
}

template<class T>
Vector4_t<T> Vector4_t<T>::cross (const Vector4_t<T>& v) const
{
    return Vector3_t<T>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x);
}

template<class T>
void Vector4_t<T>::normalize()
{
    double d = this->length();
    x /= d;
    y /= d;
    z /= d;
}

template<class T>
Vector4_t<T> Vector4_t<T>::get_normalized() const
{
    double d = this->length();
    return Vector3_t<T>(x/d, y/d, z/d, w);
}

template<class T>
void Vector4_t<T>::rotate(const double& gantry, const double& couch)
{
    double c_couch  = cos(couch);
    double s_couch  = sin(couch);
    double c_gantry = cos(gantry);
    double s_gantry = sin(gantry);

    Vector4_t<T> temp(x, y, z);
    x = temp.x*c_couch - s_couch*(temp.y*s_gantry + temp.z*c_gantry);
    y = temp.y*c_gantry - temp.z*s_gantry;
    z = temp.x*s_couch + c_couch*(temp.y*s_gantry + temp.z*c_gantry);
}

template<class T>
Vector4_t<T> Vector4_t<T>::get_rotated(const double& gantry, const double& couch) const
{
    double c_couch  = cos(couch);
    double s_couch  = sin(couch);
    double c_gantry = cos(gantry);
    double s_gantry = sin(gantry);

    Vector4_t<T> temp(x, y, z, w);
    double a_ = temp.x*c_couch - s_couch*(temp.y*s_gantry + temp.z*c_gantry);
    double b_ = temp.y*c_gantry - temp.z*s_gantry;
    double c_ = temp.x*s_couch + c_couch*(temp.y*s_gantry + temp.z*c_gantry);
    return Vector3_t<T>(a_, b_, c_, temp.w);
}

// OPERATORS AND ACCESS
template<class T>
const T& Vector4_t<T>::at(int idx) const
{
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

template<class T>
const T& Vector4_t<T>::operator [](int idx) const
{
    if (idx == 0)
        return x;
    if (idx == 1)
        return y;
    if (idx == 2)
        return z;
    if (idx == 3)
        return w;
    std::cerr << "ERROR! Index out of range!" << std::endl;
    exit(EXIT_FAILURE);
    return x;
}

template<class T> template<class U>
Vector4_t<T>& Vector4_t<T>::operator+=(const Vector4_t<U>& rhs)
{
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    w += rhs.w;
    return *this;
}
template<class T> template<class U>
Vector4_t<T>& Vector4_t<T>::operator+=(const Vector3_t<U>& rhs)
{
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
}
template<class T> template<typename U>
Vector4_t<T>& Vector4_t<T>::operator+=(const U& rhs)
{
    x += rhs;
    y += rhs;
    z += rhs;
    w += rhs;
    return *this;
}
template<class T> template<class U>
Vector4_t<T>& Vector4_t<T>::operator-=(const Vector4_t<U>& rhs)
{
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    w -= rhs.w;
    return *this;
}
template<class T> template<class U>
Vector4_t<T>& Vector4_t<T>::operator-=(const Vector3_t<U>& rhs)
{
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
}
template<class T> template<typename U>
Vector4_t<T>& Vector4_t<T>::operator-=(const U& rhs)
{
    x -= rhs;
    y -= rhs;
    z -= rhs;
    w -= rhs;
    return *this;
}
template<class T> template<class U>
Vector4_t<T>& Vector4_t<T>::operator*=(const Vector4_t<U>& rhs)
{
    x *= rhs.x;
    y *= rhs.y;
    z *= rhs.z;
    w *= rhs.w;
    return *this;
}
template<class T> template<class U>
Vector4_t<T>& Vector4_t<T>::operator*=(const Vector3_t<U>& rhs)
{
    x *= rhs.x;
    y *= rhs.y;
    z *= rhs.z;
    return *this;
}
template<class T> template<typename U>
Vector4_t<T>& Vector4_t<T>::operator*=(const U& rhs)
{
    x *= rhs;
    y *= rhs;
    z *= rhs;
    w *= rhs;
    return *this;
}
template<class T> template<class U>
Vector4_t<T>& Vector4_t<T>::operator/=(const Vector4_t<U>& rhs)
{
    x /= rhs.x;
    y /= rhs.y;
    z /= rhs.z;
    w /= rhs.w;
    return *this;
}
template<class T> template<class U>
Vector4_t<T>& Vector4_t<T>::operator/=(const Vector3_t<U>& rhs)
{
    x /= rhs.x;
    y /= rhs.y;
    z /= rhs.z;
    return *this;
}
template<class T> template<typename U>
Vector4_t<T>& Vector4_t<T>::operator/=(const U& rhs)
{
    x /= rhs;
    y /= rhs;
    z /= rhs;
    w /= rhs;
    return *this;
}

template<class T, typename U>
Vector4_t<T> operator+(Vector4_t<T> lhs, const U& rhs)
{
    return lhs += rhs;
}
template<class T, typename U>
Vector4_t<T> operator+(const U& rhs, Vector4_t<T> lhs)
{
    return lhs += rhs;
}
template<class T, typename U>
Vector4_t<T> operator+(Vector4_t<T> lhs, const Vector4_t<U>& rhs)
{
    return lhs += rhs;
}

template<class T, typename U>
Vector4_t<T> operator-(Vector4_t<T> lhs, const U& rhs)
{
    return lhs -= rhs;
}
template<class T, typename U>
Vector4_t<T> operator-(const U& rhs, Vector4_t<T> lhs)
{
    return lhs -= rhs;
}
template<class T, typename U>
Vector4_t<T> operator-(Vector4_t<T> lhs, const Vector4_t<U>& rhs)
{
    return lhs -= rhs;
}

template<class T, typename U>
Vector4_t<T> operator*(Vector4_t<T> lhs, const U& rhs)
{
    return lhs *= rhs;
}
template<class T, typename U>
Vector4_t<T> operator*(const U& rhs, Vector4_t<T> lhs)
{
    return lhs *= rhs;
}
template<class T, typename U>
Vector4_t<T> operator*(Vector4_t<T> lhs, const Vector4_t<U>& rhs)
{
    return lhs *= rhs;
}

template<class T, typename U>
Vector4_t<T> operator/(Vector4_t<T> lhs, const U& rhs)
{
    return lhs /= rhs;
}
template<class T, typename U>
Vector4_t<T> operator/(const U& rhs, Vector4_t<T> lhs)
{
    return lhs /= rhs;
}
template<class T, typename U>
Vector4_t<T> operator/(Vector4_t<T> lhs, const Vector4_t<U>& rhs)
{
    return lhs /= rhs;
}
