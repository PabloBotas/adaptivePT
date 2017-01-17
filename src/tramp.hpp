#ifndef __TRAMP_HPP__
#define __TRAMP_HPP__

#include <string>
#include <vector>

#include "spot.hpp"

class Tramp_t
{
public:
    Tramp_t();
    Tramp_t(std::string file);
    ~Tramp_t();
    void setEnergies();
    void to_file(std::string file);
    void scale(double ratio);
    std::vector<double> getWEPLs();
    void shiftEnergies(std::vector<double> ds);
    void print(unsigned int n);
    void print(unsigned int n0, unsigned int n1);

    std::vector<Spot_t> spots;
    std::vector<double> energies;

    double gigaprotons;
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
private:
    template<class T>
    T getHeaderValue(std::string line);
    template<class T>
    T getHeaderValue(std::ifstream &stream);
    std::string file;
    void read_();
};

#endif
