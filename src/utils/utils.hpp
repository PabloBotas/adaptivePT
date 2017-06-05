#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <string>
#include <vector>
#include "special_types.hpp"

namespace utils
{
    std::string toLower(std::string s);
    bool ends_with_string(std::string const& str, std::string const& what);
    std::string get_file_extension(std::string const& str);
    std::string replace_substring(std::string const& str,
                              std::string const& what,
                              std::string const& to_replace);
    std::string run_command(const std::string cmd);
    std::string to_location_str(const Vector3_t<float>& p, const bool b);
    void flip_positions_X(std::vector< Vector4_t<float> >& vec,
                          const CT_Dims_t dims);
    void run_plastimatch_probe(const std::vector< Vector4_t<float> >& p,
                                              const std::string vf);
    void run_plastimatch_probe(const std::vector< Vector3_t<float> >& p,
                                              const std::string vf);
    void cm_to_mm(std::vector< Vector4_t<float> >& v);
}

#endif