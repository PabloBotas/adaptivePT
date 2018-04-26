#include "utils.hpp"
#include "volume.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <stdlib.h>
#include <fstream>
#include <array>
#include <regex>
#include <cerrno>
#include <omp.h>

template <class T>
std::vector<T> utils::join_vectors (const std::vector<T>& a, const std::vector<T>& b)
{
    std::vector<T> out(a.size() + b.size());
    std::copy(a.begin(), a.end(), out.begin());
    std::copy(b.begin(), b.end(), out.begin()+a.size());
    return out;
}
template std::vector<std::string> utils::join_vectors<std::string>(const std::vector<std::string>&,
                                                                   const std::vector<std::string>&);


template <class T>
std::vector<T> utils::join_vectors (const std::vector<T>& a, const std::vector<T>& b,
                                    const std::vector<T>& c)
{
    std::vector<T> out(a.size() + b.size() + c.size());
    std::copy(a.begin(), a.end(), out.begin());
    std::copy(b.begin(), b.end(), out.begin()+a.size());
    std::copy(c.begin(), c.end(), out.begin()+a.size()+b.size());
    return out;
}
template std::vector<std::string> utils::join_vectors<std::string>(const std::vector<std::string>&,
                                                                   const std::vector<std::string>&,
                                                                   const std::vector<std::string>&);


template <class T>
std::vector<T> utils::or_vectors (const std::vector<T>& a, const std::vector<T>& b)
{
    if (a.size() != b.size()) {
        std::cerr << "ERROR! Vector sizes passed to or_vectors are not consistent!" << std::endl;
        std::cerr << a.size() << " != " << b.size() << std::endl;
        exit(EXIT_FAILURE); 
    }
    std::vector<T> out = a;
    for (size_t i = 0; i < a.size(); ++i) {
        if (out[i] > 0.5)
            continue;
        out[i] = b[i];
    }
    return out;
}
template std::vector<float> utils::or_vectors<float>(const std::vector<float>&,
                                                     const std::vector<float>&);


template <class T>
std::vector<T> utils::or_vectors (const std::vector<T>& a, const std::vector<T>& b,
                                  const std::vector<T>& c)
{
    if ((a.size() != b.size()) || (a.size() != c.size()) || (b.size() != c.size()))  {
        std::cerr << "ERROR! Vector sizes passed to or_vectors are not consistent!" << std::endl;
        std::cerr << a.size() << ", " << b.size() << ", " << c.size() << std::endl;
        exit(EXIT_FAILURE); 
    }
    std::vector<T> out = a;
    for (size_t i = 0; i < a.size(); ++i) {
        if (out[i] > 0.5)
            continue;
        out[i] = b[i] + c[i] > 0.5 ? 1 : 0;
    }    
    return out;
}
template std::vector<float> utils::or_vectors<float>(const std::vector<float>&,
                                                     const std::vector<float>&,
                                                     const std::vector<float>&);


void utils::check_fs(const std::ofstream& fs, std::string f, std::string msg)
{
    if (fs.fail()) {
        std::cerr << "Can not open file \"" << f << "\" " << msg << std::endl;
        exit(EXIT_FAILURE);
    }
}
void utils::check_fs(const std::ifstream& fs, std::string f, std::string msg)
{
    if (fs.fail()) {
        std::cerr << "Can not open file \"" << f << "\" " << msg << std::endl;
        exit(EXIT_FAILURE);
    }
}

template <class T>
int utils::sgn(T& v)
{
    return (v > T(0)) - (v < T(0));
}
template int utils::sgn<float>(float&);

template <class T>
int utils::sgn_plus(T& v)
{
    return (v >= T(0)) - (v < T(0));
}
template int utils::sgn_plus<float>(float&);

template <class T>
int utils::sgn_minus(T& v)
{
    return (v > T(0)) - (v <= T(0));
}
template int utils::sgn_minus<float>(float&);

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

std::string utils::get_parent_path(const std::string& str)
{
    size_t pos = str.rfind('/');
    if (pos != std::string::npos) {
        return str.substr(0, str.rfind('/'));
    } else {
        return std::string(".");
    }
}

std::string utils::get_full_path(const std::string& str)
{
    char *actualpath = realpath(str.c_str(), NULL);
    std::string path(actualpath);
    free(actualpath);
    return path;
}

std::string utils::get_file_name(const std::string& str)
{
    size_t pos = str.rfind('/');
    if (pos != std::string::npos) {
        return str.substr(pos+1);
    } else {
        return str;
    }
}

std::string utils::get_full_parent_path(const std::string& str)
{
    return get_parent_path(get_full_path(str));
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
                                  std::string const& new_substr)
{
    std::string out = str;
    replace_string_inplace(out, substr, new_substr);
    return out;
}

void utils::replace_string_inplace(std::string& str,
                                   const std::string& substr,
                                   const std::string& new_substr)
{
    size_t pos = str.find(substr);
    // If the new substring contains the substring, break after one replacement
    bool not_recursive = new_substr.find(substr) != std::string::npos ? true : false;
    while (pos != std::string::npos) {
        str.replace(pos, substr.length(), new_substr);
        pos = str.find(substr, pos);
        if (not_recursive)
            break;
    }
}


void utils::copy_file(const std::string& in,
                      const std::string& out)
{
    std::ifstream src(in);
    if (!src) {
        std::cerr << "Can't open \"" + in + "\" file to copy" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::ofstream dst(out);
    if (!dst) {
        std::cerr << "Can't open \"" + out + "\" file to copy" << std::endl;
        exit(EXIT_FAILURE);
    }
    dst << src.rdbuf();
}


void utils::copy_replace_in_file(const std::string& in,
                                 const std::string& out,
                                 std::map<std::string, std::string> mymap)
{
    std::ifstream src(in);
    if (!src) {
        std::cerr << "Can't open \"" + in + "\" file to copy" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::ofstream dst(out);
    if (!dst) {
        std::cerr << "Can't open \"" + out + "\" file to copy" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string text;
    for (char ch; src.get(ch); text.push_back(ch)) {
    }
    for(auto it = mymap.begin(); it != mymap.end(); ++it) {
        replace_string_inplace(text, it->first, it->second);
    }
    dst << text;
}


Volume_t utils::read_masks (const std::vector<std::string>& v, const float threshold,
                            const std::vector<int> mask_importances) {
    // Read structure volume
    Volume_t vol(v.front(), Volume_t::Source_type::MHA);
    for (size_t f=1; f < v.size(); ++f) {
        Volume_t temp(v.at(f), Volume_t::Source_type::MHA);
        if (temp.nElements != vol.nElements ||
            temp.n.x != vol.n.x || temp.n.y != vol.n.y || temp.n.z != vol.n.z ||
            temp.d.x != vol.d.x || temp.d.y != vol.d.y || temp.d.z != vol.d.z) {
            std::cerr << "ERROR! Dimensions of passed mask " << v.at(f);
            std::cerr << " differ from first passed mask. Plastimatch says: " << std::endl;
            std::string cmd("plastimatch header " + v.at(0) + " 2>&1");
            std::string temp = utils::run_command(cmd);
            std::cerr << v.front() << ":" << std::endl;
            std::cerr << temp << std::endl;
            cmd = "plastimatch header " + v.at(f) + " 2>&1";
            temp = utils::run_command(cmd);
            std::cerr << v.at(f) << ":" << std::endl;
            std::cerr << temp << std::endl;
            exit(EXIT_FAILURE);
        }
        // The masks will be passed as int to the device. To prevent unwanted flooring due to
        // precision, I add 0.5.
        float val = mask_importances.size() ? mask_importances.at(f) : 2*threshold;
        #pragma omp parallel for
        for (uint i = 0; i < temp.nElements; ++i) {
            if (temp.data.at(i) > threshold && vol.data.at(i) <= threshold) {
                vol.data.at(i) = val;
            }
        }
    }

    return vol;
}


void utils::append_to_file(const std::string& f,
                           std::string append)
{
    std::ofstream src(f, std::ios::ate | std::ios::app);
    if (!src) {
        std::cerr << "Can't open " + f + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
    src << append;
}


std::string utils::run_command(const std::string cmd, int* return_code)
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
    auto deleter = [return_code] (std::FILE* p) {
        if (return_code)
            *return_code = pclose(p)/256;
        else
            pclose(p);
    };
    std::shared_ptr<std::FILE> pipe(popen(cmd.c_str(), "r"), deleter);
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
utils::subset_vector< Vector4_t<float>>(std::vector< Vector4_t<float>>&,
                                        const std::vector< Vector4_t<float>>&,
                                        const size_t, const size_t);

double utils::range_from_energy(const double& e, bool conv)
{
    // energy in MeV range in cm of water
    double E = e;
    if (conv)
        E /= 1e6;
    double r_e[6] = {1.53073818e-12, -4.52094359e-10, -5.83436734e-07, 6.95245227e-04,
                     1.61622821e-02, -3.55845392e-01};
    double R = r_e[5] + r_e[4]*E + r_e[3]*std::pow(E, 2) + r_e[2]*std::pow(E, 3) +
               r_e[1]*std::pow(E, 4) + r_e[0]*std::pow(E, 5);
    return R;
}

double utils::energy_from_range(const double& r)
{
    // energy in MeV range in cm of water
    double e_r[8] = {6.10690676e-09, -1.02370849e-06, 7.18732504e-05, -2.76315988e-03,
                     6.43961648e-02, -9.85218453e-01, 1.49352594e+01, 2.27490185e+01};
    double E = e_r[7] + e_r[6]*r + e_r[5]*std::pow(r, 2) + e_r[4]*std::pow(r, 3) +
               e_r[3]*std::pow(r, 4) + e_r[2]*std::pow(r, 5) + e_r[1]*std::pow(r, 6) +
               e_r[0]*std::pow(r, 7);
    // std::cerr << "energy_from_range " << r << " " << E << std::endl;
    return E*1e6;
}

