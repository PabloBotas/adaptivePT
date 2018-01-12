#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include "special_types.hpp"
#include "vector4.hpp"
#include "volume.hpp"

#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace utils
{
    std::string toLower(std::string s);
    bool starts_with_string(std::string const& str, std::string const& what);
    bool ends_with_string(std::string const& str, std::string const& what);
    std::string get_full_path(const std::string& str);
    std::string get_full_parent_path(const std::string& str);
    std::string get_parent_path(const std::string& str);
    std::string get_file_name(std::string const& path);
    std::string get_file_extension(std::string const& str);
    std::string remove_file_extension(std::string const& str);
    std::string replace_string(std::string const& str,
                               std::string const& substr,
                               std::string const& to_replace);
    void replace_string_inplace(std::string& str,
                                std::string const& substr,
                                std::string const& to_replace);
    void copy_file(const std::string& in,
                   const std::string& out);
    void copy_replace_in_file(const std::string& in,
                              const std::string& out,
                              std::map<std::string, std::string> mymap);
    void append_to_file(const std::string& f,
                        std::string append);
    Volume_t read_masks (const std::vector<std::string>& v, const float threshold = 0.5);
    std::string run_command(const std::string cmd);
    void cm_to_mm(Array4<float>& v);
    Vector3_t<float> closest_point(const Vector3_t<float>& vec,
                                   const Vector3_t<float>& vec_p,
                                   const Vector3_t<float>& p);
    void check_fs(const std::ofstream& ofs, std::string f, std::string msg);
    void check_fs(const std::ifstream& ofs, std::string f, std::string msg);
    template<class T>
    void subset_vector(std::vector<T>& vout, const std::vector<T>& v,
                       const size_t offset_a, const size_t offset_b);
}

#endif