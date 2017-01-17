#include "data_table.hpp"

#include "spline.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>


//----------------------------------------------------
// Constructors --------------------------------------
//----------------------------------------------------
DataTable_t::DataTable_t()
{
}

DataTable_t::~DataTable_t()
{
}

DataTable_t::DataTable_t(std::string colname, std::vector<double> vec)
{
    append(colname, vec);
}

DataTable_t::DataTable_t(std::string colname1, std::vector<double> x,
                         std::string colname2, std::vector<double> y)
{
    append(colname1, x);
    append(colname2, y);
}

DataTable_t::DataTable_t(size_t n, size_t m)
{
    resize(Dim_t::Y, m);
    resize(Dim_t::X, n);
}

//----------------------------------------------------
// Constructors --------------------------------------
//----------------------------------------------------
void DataTable_t::append(std::string colname, std::vector<double> vec)
{
    cols.push_back(colname);
    data.push_back(vec);
}

void DataTable_t::resize(Dim_t dim, size_t sz)
{
    switch(dim)
    {
        case Dim_t::X:
            if (data.empty())
            {
                data.resize(1);
            }
            for (size_t i = 0; i < data.size(); i++)
            {
                data.at(i).resize(sz);
            }
            break;
        case Dim_t::Y:
            data.resize(sz);
            break;
    }
}

//----------------------------------------------------
// Acess ---------------------------------------------
//----------------------------------------------------
std::vector<double>& DataTable_t::at(size_t i)
{
    return data.at(i);
}

double DataTable_t::at(size_t i, size_t j)
{
    return data.at(i).at(j);
}

double DataTable_t::at(std::string colname, size_t j)
{
    size_t colnum = getColNumber(colname);

    return data.at(colnum).at(j);
}


//----------------------------------------------------
// Queries -------------------------------------------
//----------------------------------------------------
size_t DataTable_t::getColNumber(std::string colname)
{
    std::vector<std::string>::iterator it = std::find(cols.begin(), cols.end(), colname);
    if (std::distance(it, cols.end()) == 0)
    {
        std::cerr << "ERROR! Column " << colname << " not found in data table." << std::endl;
        std::copy(cols.begin(), cols.end(), std::ostream_iterator<std::string>(std::cerr, " "));
    }
    size_t colnum = std::distance(cols.begin(), it);

    return colnum;
}

size_t DataTable_t::getNumberRows()
{
    return data.size() > 0 ? data.front().size() : 0;
}

size_t DataTable_t::getNumberCols()
{
    return data.size() > 0 ? data.front().size() : 0;
}

double DataTable_t::getVal(std::string col_out_name, std::string col_in_name, double x)
{
    size_t col_out = getColNumber(col_out_name);
    size_t col_in  = getColNumber(col_in_name);
    std::vector<double>::iterator begin = data.at(col_in).begin();
    std::vector<double>::iterator end   = data.at(col_in).end();
    std::vector<double>::iterator it = std::upper_bound(begin, end, x);
    size_t idx = std::distance(begin, it);
    double val = *it;
    
    double out;
    if( round(100.*val) != round(100.*x) )
    {
        // Value not exact in x, call spline
        Spline_t spline;
        spline.set_points(data.at(col_in), data.at(col_out));
        out = spline(x);
    }
    else
    {
        out = data.at(col_out).at(idx);
    }

    return out;
}



