#ifndef __SPOT_HPP__
#define __SPOT_HPP__

#include <string>
#include <vector>

struct Spot_t
{
    Spot_t();
    Spot_t(std::string line);
    Spot_t(float e, float w, float x, float y);
    ~Spot_t();
    float x;
    float y;
    float e;
    float w;

    void ShiftEnergy(float d);
    friend std::ostream& operator<<(std::ostream& os, const Spot_t& s);
};

#endif
