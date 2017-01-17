#ifndef __SPOT_HPP__
#define __SPOT_HPP__

#include <string>
#include <vector>

struct Spot_t
{
    Spot_t();
    Spot_t(std::string line);
    Spot_t(double e, double w, double x, double y);
    ~Spot_t();
    std::string x;
    std::string y;
    double e;
    double w;

    void ShiftEnergy(double d);
    friend std::ostream& operator<<(std::ostream& os, const Spot_t& s);
};

#endif
