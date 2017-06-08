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
    std::cout << "Running command:";

    if(cmd.size() > 256)
    {
        std::cout << " (trimmed)" << std::endl;
        std::cout.write(&cmd[0], 256);
        std::cout << " ...//... ";
        size_t size = cmd.size();
        std::cout.write(&cmd[size-128], 128) << std::endl;
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


void 
utils::run_plastimatch_probe(std::vector< Vector4_t<float> >& p,
                             std::string vf_file)
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
    if ( ! (utils::ends_with_string(vf_file, ".mha") || 
            utils::ends_with_string(vf_file, ".mhd")) )
    {
        std::string ext = utils::get_file_extension(vf_file);
        std::string trans_cmd;
        trans_cmd = "plastimatch xf-convert --input " + vf_file;
        vf_file = utils::replace_substring(vf_file, ext, "mha");
        trans_cmd += " --output " + vf_file +
                     " --output-type vf";

        std::string temp = utils::run_command(trans_cmd);
        std::cout << temp << std::endl;
    }

    std::string cmd = "plastimatch probe --location \"";
    cmd += locations + "\" ";
    cmd += vf_file;

    // Run command
    std::string stdout = utils::run_command(cmd);
    // std::cout << stdout << std::endl;

    std::vector< Vector3_t<float> > vf = utils::get_vf_from_stdout(stdout);

    for (size_t i = 0; i < vf.size(); i++)
    {
        p.at(i).x += vf.at(i).x;
        p.at(i).y += vf.at(i).y;
        p.at(i).z += vf.at(i).z;
    }
}

std::vector< Vector3_t<float> >
utils::get_vf_from_stdout(std::string str)
{
    std::vector< Vector3_t<float> > v;

    std::stringstream ss_per_line(str);
    std::string line;

    if (!str.empty())
    {
        while(std::getline(ss_per_line, line, '\n'))
        {
            std::istringstream ss_per_space(line);

            Vector3_t<float> this_vf;
            while (ss_per_space)
            {
                std::string dummy;
                ss_per_space >> dummy;
                ss_per_space >> dummy >> dummy >> dummy;
                ss_per_space >> dummy >> dummy >> dummy;
                ss_per_space >> this_vf.x >> this_vf.y >> this_vf.z;
            }
            v.push_back(Vector3_t<float>(this_vf.x, this_vf.y, this_vf.z));
        }
    }

    return v;
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
