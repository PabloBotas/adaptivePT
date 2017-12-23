#include "utils.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <array>
#include <regex>
#include <cerrno>

void utils::check_fs(const std::ofstream& fs, std::string f, std::string msg)
{
    if ( !fs.is_open() ) {
        std::cerr << "Can not open file " << f << " " << msg << std::endl;
        exit(EXIT_FAILURE);
    }
}
void utils::check_fs(const std::ifstream& fs, std::string f, std::string msg)
{
    if ( !fs.is_open() ) {
        std::cerr << "Can not open file " << f << " " << msg << std::endl;
        exit(EXIT_FAILURE);
    }
}


std::string utils::toLower(std::string s)
{
    transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

bool utils::starts_with_string(std::string const& str,
                               std::string const& what)
{
    return str.size() >= what.size() &&
        equal(what.begin(), what.end(), str.begin());
}

bool utils::ends_with_string(std::string const& str,
                             std::string const& what)
{
    return what.size() <= str.size() &&
           str.find(what, str.size() - what.size()) != str.npos;
}

std::string utils::get_parent_directory(const std::string& str)
{
    return str.substr(0, str.rfind('/'));
}

std::string utils::get_file_extension(std::string const& str)
{
    return str.substr(str.rfind('.') + 1);
}

std::string utils::remove_file_extension(std::string const& str)
{
    return str.substr(0, str.rfind('.'));
}

std::string utils::replace_string(std::string const& str,
                                  std::string const& substr,
                                  std::string const& to_replace)
{
    std::string out = str;
    replace_string(out, substr, to_replace);
    return out;
}

void utils::replace_string_inplace(std::string& str,
                                   const std::string& substr,
                                   const std::string& to_replace)
{
    size_t pos = str.find(substr);
    while (pos != std::string::npos) {
        str.replace(pos, substr.length(), to_replace);
        pos = str.find(substr, pos);
    }
}

void utils::copy_file(const std::string& in,
                      const std::string& out,
                      std::map<std::string, std::string> mymap)
{
    std::ifstream src(in);
    if (!src) {
        std::cerr << "Can't open " + in + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::ofstream dst(out);
    if (!dst) {
        std::cerr << "Can't open " + out + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!mymap.empty()) {
        std::string text;
        for (char ch; src.get(ch); text.push_back(ch)) {
        }
        for(std::map<std::string, std::string>::iterator it = mymap.begin(); it != mymap.end(); ++it) {
            replace_string(text, it->first, it->second);
        }
        dst << text;
    } else {
        dst << src.rdbuf();
    }
}


std::string utils::run_command(const std::string cmd)
{
    std::cout << std::endl << "Running command";

    std::string cmd_trimmed;
    const std::string* to_output = &cmd;
    if (cmd.size() > 160) {
        std::cout << " (trimmed)";
        cmd_trimmed = cmd.substr(0, 50) + " ...//... " + cmd.substr(cmd.size()-50, 50);
        to_output = &cmd_trimmed;
    }

    std::cout << ": " << std::endl << *to_output << std::endl;


    std::array<char, 512> buffer;
    std::string stdout;
    std::shared_ptr<std::FILE> pipe(popen(cmd.c_str(), "r"), pclose);
    if (!pipe) {
        std::cerr << std::endl << "Could not launch command:" << std::endl;
        std::cerr << *to_output << std::endl;
        std::cerr << "Error number: " << errno << std::endl;
        std::perror("popen() failed");
        exit(EXIT_FAILURE);
    }
    while (!feof(pipe.get())) {
        if (std::fgets(buffer.data(), 512, pipe.get()) != NULL)
            stdout += buffer.data();
    }
    return stdout;
}


void utils::cm_to_mm(Array4<float>& vec)
{
    const int CM2MM = 10;
    for (size_t i = 0; i < vec.size(); i++) {
        vec[i].x *= CM2MM;
        vec[i].y *= CM2MM;
        vec[i].z *= CM2MM;
    }
}

Vector3_t<float> utils::closest_point(const Vector3_t<float>& vec,
                                      const Vector3_t<float>& vec_p,
                                      const Vector3_t<float>& p)
{
    // Closest point along line defined by vec_p and vec from point p
    return vec_p + (vec.dot(p-vec_p)/vec.length2())*vec;
}


template<class T>
void utils::subset_vector(std::vector<T>& vout, const std::vector<T>& v,
                          const size_t offset_a, const size_t offset_b)
{
    typename std::vector<T>::const_iterator a = v.begin() + offset_a;
    typename std::vector<T>::const_iterator b = v.begin() + offset_b;
    vout.assign(a, b);
}
template void
utils::subset_vector<float>(std::vector<float>&, const std::vector<float>&,
                             const size_t, const size_t);
template void
utils::subset_vector< Vector4_t<float> >(std::vector< Vector4_t<float> >&,
                                          const std::vector< Vector4_t<float> >&,
                                          const size_t, const size_t);
