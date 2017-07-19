#include "spot.hpp"

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

Spot_t::~Spot_t()
{
}

Spot_t::Spot_t()
{
}

Spot_t::Spot_t(double e_, double w_, double x_, double y_)
{
    e = e_;
    w = w_;
    x = x_;
    y = y_;
}

Spot_t::Spot_t(std::string line)
{
    e = 0;
    w = 0;
    x = 0;
    y = 0;
    std::istringstream ss(line);
    try {
        ss >> e >> x >> y >> w;
    }
    catch (const std::exception &exc)
    {
        std::cerr << "ERROR! while parsing spot: \n";
        std::cerr << line << std::endl;
        std::cerr << exc.what();
        exit(EXIT_FAILURE);
    }
}


std::ostream& operator<<(std::ostream& os, const Spot_t& s)
{
    std::string delim_e = s.e < 100.f ? " " : "";
    std::string delim_x = " ";
    if (s.x > 0)
        delim_x += " ";
    if (std::abs(s.x) < 100.f)
        delim_x += " ";
    if (std::abs(s.x) < 10.f)
        delim_x += " ";
    std::string delim_y = " ";
    if (s.y > 0)
        delim_y += " ";
    if (std::abs(s.y) < 100.f)
        delim_y += " ";
    if (std::abs(s.y) < 10.f)
        delim_y += " ";
    std::string delim_w = "   ";
    os << delim_e << s.e << delim_x << s.x << delim_y << s.y << delim_w << s.w;
    return os;  
}

