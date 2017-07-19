#include "utils.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
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

    std::array<char, 512> buffer;
    std::string stdout;
    std::shared_ptr<std::FILE> pipe(popen(cmd.c_str(), "r"), pclose);
    if (!pipe)
    {
        std::cerr << "ERROR! popen() failed" << std::endl;
        exit(EXIT_FAILURE);
    }
    while (!feof(pipe.get())) {
        if (std::fgets(buffer.data(), 512, pipe.get()) != NULL)
            stdout += buffer.data();
    }
    return stdout;
}


void utils::cm_to_mm(Array4<double>& vec)
{
    const int CM2MM = 10;
    for (size_t i = 0; i < vec.size(); i++)
    {
        vec[i].x *= CM2MM;
        vec[i].y *= CM2MM;
        vec[i].z *= CM2MM;
    }
}

Vector3_t<double> utils::intersect(const Vector3_t<double>& a,
                                   const Vector3_t<double>& u,
                                   const Vector3_t<double>& b,
                                   const Vector3_t<double>& v)
{
    // Check they are coplanar
    Vector3_t<double> cross = u.cross(v);
    if (cross.length() == 0.f)
    {
        std::cout << "WARNING! Intersect received non-coplanar lines!" << std::endl;
        std::cout << "Detected: F: " << __FILE__ << ". L: " << __LINE__ << std::endl;
        return Vector3_t<double>(NAN, NAN, NAN);
    }

    // Some direction components may be zero
    // I find the largest one to be used as denominator
    // Then I just select the next one in the container
    // (operator[] has recursive nature in VectorN_t class)
    int ind_max = (u.x > u.y) ? (u.x > u.z ? 0 : 2) : (u.y > u.z ? 1 : 2);
    int ind_min = (u.x < u.y) ? (u.x < u.z ? 0 : 2) : (u.y < u.z ? 1 : 2);
    double u_ratio = u[ind_min]/u[ind_max];
    
    double s = ((a[ind_min]-b[ind_min]) + (b[ind_max]-a[ind_max])*u_ratio) / (v[ind_min] - v[ind_max]*u_ratio); // line parameter of b,v
    double t = ((b[ind_max]-a[ind_max]) + s*v[ind_max])/u[ind_max];
    Vector3_t<double> A = a + t*u;
    Vector3_t<double> B = b + s*v;

    std::cout << "ind_max: " << ind_max << std::endl;
    std::cout << "ind_min: " << ind_min << std::endl;
    std::cout << "a: " << a.x << " " << a.y << " " << a.z << std::endl;
    std::cout << "u: " << u.x << " " << u.y << " " << u.z << std::endl;
    std::cout << "b: " << b.x << " " << b.y << " " << b.z << std::endl;
    std::cout << "v: " << v.x << " " << v.y << " " << v.z << std::endl;
    std::cout << "t: " << t << std::endl;
    std::cout << "s: " << s << std::endl;
    std::cout << "A: " << A.x << " " << A.y << " " << A.z << std::endl;
    std::cout << "B: " << B.x << " " << B.y << " " << B.z << std::endl;

    return A;
}

