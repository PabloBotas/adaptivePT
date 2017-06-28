#ifndef __TRAMP_HPP__
#define __TRAMP_HPP__

#include <string>
#include <vector>

#include "spot.hpp"
#include "vector4.hpp"

class Tramp_t
{
public:
    Tramp_t();
    Tramp_t(std::string file);
    Tramp_t(std::string file, std::string machine);
    ~Tramp_t();
    void setEnergies();
    void to_file(std::string file, std::string dir = std::string());
    void scale(float ratio);
    void setWEPLs();
    std::vector<float> getWEPLs();
    void energy_to_internal();
    void shift_energies(const std::vector<float>& e_, bool units = true);
    void set_pos(const std::vector< Vector4_t<float> > p);
    void print(unsigned int n);
    void print(unsigned int n0, unsigned int n1);
    void defaults();
    void read_file_header(std::string f);

    std::vector<Spot_t> spots;
    std::vector<float> energies;
    std::vector<float> energies_internal;
    std::vector<float> wepls;

    std::string machine;

    float gigaprotons;
    unsigned int nspots;
    std::string patient_id;
    std::string patient_first_name;
    std::string patient_middle_initial;
    std::string patient_last_name;
    std::string astroid_id;
    std::string course_name;
    std::string beam_name;
    std::string gantry;
    std::string couch_rotation;

    float z;
    float zeff;
    
private:
    bool internal_energy_set;
    template<class T>
    T getHeaderValue(std::string line);
    template<class T>
    T getHeaderValue(std::ifstream &stream);
    std::string file;
    void read_();
};

std::string toLower(std::string s);


#endif
