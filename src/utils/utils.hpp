#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include "special_types.hpp"

#include <string>
#include <vector>

namespace utils
{
    std::string toLower(std::string s);
    bool ends_with_string(std::string const& str, std::string const& what);
    std::string get_file_extension(std::string const& str);
    std::string replace_substring(std::string const& str,
                              std::string const& what,
                              std::string const& to_replace);
    std::string run_command(const std::string cmd);
    void cm_to_mm(Array4<float>& v);
}

#endif