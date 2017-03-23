#include "utils.hpp"

#include <string>
#include <algorithm>

std::string utils::toLower(std::string s)
{
    transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}
