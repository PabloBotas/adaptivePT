#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include "special_types.hpp"
#include "vector4.hpp"

#include <fstream>
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
    Vector4_t<float> rotate(const Vector4_t<float>& p, const float& gantry, const float& couch);
    void check_fs(const std::ofstream& ofs, std::string f, std::string msg);
    void check_fs(const std::ifstream& ofs, std::string f, std::string msg);
}

#endif