#include "special_types.hpp"

#include <vector>

template<class T>
Vector_t<T>::Vector_t()
{
    x = 0;
    y = 0;
    z = 0;
}

template<class T>
Vector_t<T>::Vector_t(const std::vector<T> &v)
{
    x = v.at(0);
    y = v.at(1);
    z = v.at(2);
}

template<class T>
Vector_t<T>::Vector_t(const Vector_t<T> &obj)
{
    x = obj.x;
    y = obj.y;
    z = obj.z;
}

template class Vector_t<int>;
template class Vector_t<unsigned int>;
template class Vector_t<float>;
template class Vector_t<double>;
