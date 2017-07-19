#ifndef __ENERGY_RANGE_HPP__
#define __ENERGY_RANGE_HPP__

#include "data_table.hpp"

#include <string>
#include <vector>

class Energy_Range_Calculator_t
{
public:
    enum Direction_t {FromEtoR, FromRtoE};
    Direction_t dir;
    Energy_Range_Calculator_t(Direction_t dir = FromEtoR);
    double calculate(double x, Direction_t dir = FromEtoR);
    double operator()(double& x);
private:
    DataTable_t table;
};

#endif
