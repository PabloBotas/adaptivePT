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
    void cm_to_mm(Array4<double>& v);
    Vector3_t<double> closest_point(const Vector3_t<double>& vec,
                                    const Vector3_t<double>& vec_p,
                                    const Vector3_t<double>& p);
    void check_fs(const std::ofstream& ofs, std::string f, std::string msg);
    void check_fs(const std::ifstream& ofs, std::string f, std::string msg);
    template<class T>
    void subset_vector(std::vector<T>& vout, const std::vector<T>& v,
                       const size_t offset_a, const size_t offset_b);
}

#endif