#ifndef __ENERGY_RANGE_HPP__
#define __ENERGY_RANGE_HPP__

#include "data_table.hpp"

#include <string>
#include <vector>

class Energy_Range_Calculator_t
{
public:
    enum Dir_t {FromEtoR, FromRtoE};
    Dir_t dir;
    Energy_Range_Calculator_t(Dir_t dir=FromEtoR);
    double calculate(double x, Dir_t dir=FromEtoR);
    double operator()(double& x);
private:
    DataTable_t table;
};

#endif