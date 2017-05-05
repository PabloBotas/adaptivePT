#include "spot.hpp"

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

Spot_t::~Spot_t()
{
}

Spot_t::Spot_t(float e_, float w_, float x_, float y_)
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
    os << s.e << " " << s.x << " " << s.y << " " << s.w;
    return os;  
}

