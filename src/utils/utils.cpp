#include "utils.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <array>
#include <regex>

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
    std::cerr << "Running command:" << std::endl;
    std::cerr << cmd << std::endl;

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


std::string utils::to_location_str(const Vector3_t<float>& p, const bool last)
{
    std::string s;
    s += " " + std::to_string(p.x) +
         " " + std::to_string(p.y) +
         " " + std::to_string(p.z);
    if (!last)
        s += ";";
    return s;
}


void utils::run_plastimatch_probe(const std::vector< Vector3_t<float> >& p,
                                  std::string vf)
{
   //  plastimatch probe --location "0 0 0; 0.5 0.5 0.5; 1 1 1" infile.nrrd

    // Build command
    std::string locations;
    size_t elements = p.size();
    for (size_t i = 0; i < elements; i++)
    {
        locations += utils::to_location_str(p.at(i), i+1 == elements);
    }

    // Get vector field image
    if ( ! (utils::ends_with_string(vf, ".mha") || 
            utils::ends_with_string(vf, ".mhd")) )
    {
        std::string ext = utils::get_file_extension(vf);
        std::string trans_cmd;
        trans_cmd = "plastimatch xf-convert --input " + vf;
        vf = utils::replace_substring(vf, ext, "mha");
        trans_cmd += " --output " + vf +
                     " --output-type vf";

        std::string temp = utils::run_command(trans_cmd);
        std::cout << temp << std::endl;
    }

    std::string cmd = "plastimatch probe --location \"";
    cmd += locations + "\" ";
    cmd += vf;

    // Run command
    std::string stdout = utils::run_command(cmd);
    std::cout << stdout << std::endl;
}

void utils::run_plastimatch_probe(const std::vector< Vector4_t<float> >& p,
                                  std::string vf)
{
   std::vector< Vector3_t<float> > temp(p.size());
   for (size_t i = 0; i < p.size(); i++)
   {
       temp.at(i) = p.at(i);
   }
   utils::run_plastimatch_probe(temp, vf);
}


void utils::flip_positions_X(std::vector< Vector4_t<float> >& vec,
                             const CT_Dims_t dims)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        vec.at(i).x = dims.n.x*dims.d.x - vec.at(i).x;
    }
}

void utils::cm_to_mm(std::vector< Vector4_t<float> >& vec)
{
    const int CM2MM = 10;
    for (size_t i = 0; i < vec.size(); i++)
    {
        vec[i].x *= CM2MM;
        vec[i].y *= CM2MM;
        vec[i].z *= CM2MM;
    }
}
