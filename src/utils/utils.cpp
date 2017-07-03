#include "utils.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <array>
#include <regex>

void utils::check_fs(const std::ofstream& fs, std::string f, std::string msg)
{
    if( !fs.is_open() )
    {
        std::cerr << "Can not open file " << f << " " << msg << std::endl;
        exit(EXIT_FAILURE);
    }
}
void utils::check_fs(const std::ifstream& fs, std::string f, std::string msg)
{
    if( !fs.is_open() )
    {
        std::cerr << "Can not open file " << f << " " << msg << std::endl;
        exit(EXIT_FAILURE);
    }
}


std::string utils::toLower(std::string s)
{
    transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

bool utils::ends_with_string(std::string const& str,
                      std::string const& what)
{
    return what.size() <= str.size() &&
           str.find(what, str.size() - what.size()) != str.npos;
}

std::string utils::get_file_extension(std::string const& str)
{
    return str.substr(str.rfind('.') + 1);
}

std::string utils::replace_substring(std::string const& str,
                                     std::string const& what,
                                     std::string const& to_replace)
{
    std::string out = str;
    for (size_t pos = 0; ; pos += to_replace.length())
    {
        // Locate the substring to replace
        pos = str.find(what, pos );
        if( pos == std::string::npos ) break;
        // Replace by erasing and inserting
        out.erase(pos, what.length());
        out.insert( pos, to_replace );
    }
    return out;
}

std::string utils::run_command(const std::string cmd)
{
    std::cout << "Running command:";

    if(cmd.size() > 160)
    {
        std::cout << " (trimmed)" << std::endl;
        std::cout.write(&cmd[0], 75);
        std::cout << " ...//... ";
        size_t size = cmd.size();
        std::cout.write(&cmd[size-75], 75) << std::endl;
    }
    else
    {
        std::cout << std::endl << cmd << std::endl;
    }

    int const buff_length = 128;
    std::array<char, buff_length> buffer;
    std::string stdout;
    std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
    if (!pipe)
    {
        std::cerr << "ERROR! popen() failed" << std::endl;
        exit(EXIT_FAILURE);
    }
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), buff_length, pipe.get()) != NULL)
            stdout += buffer.data();
    }
    return stdout;
}


void utils::cm_to_mm(Array4<float>& vec)
{
    const int CM2MM = 10;
    for (size_t i = 0; i < vec.size(); i++)
    {
        vec[i].x *= CM2MM;
        vec[i].y *= CM2MM;
        vec[i].z *= CM2MM;
    }
}

Vector4_t<float> utils::rotate(const Vector4_t<float>& p, const float& gantry, const float& couch)
{
    float c_couch  = cos(couch);
    float s_couch  = sin(couch);
    float c_gantry = cos(gantry);
    float s_gantry = sin(gantry);

    Vector4_t<float> res;
    res.x = p.x*c_couch - s_couch*(p.y*s_gantry + p.z*c_gantry);
    res.y = p.y*c_gantry - p.z*s_gantry;
    res.z = p.x*s_couch + c_couch*(p.y*s_gantry + p.z*c_gantry);
    res.w = p.w;

    return res;
}
