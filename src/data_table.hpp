#ifndef __DATA_TABLE_HPP__
#define __DATA_TABLE_HPP__

#include <string>
#include <vector>

class DataTable_t
{
public:
    enum Dim_t { X=1, Y=2 };
    // Constructors
    DataTable_t();
    ~DataTable_t();
    DataTable_t(std::string colname, std::vector<double> x);
    DataTable_t(std::string colname1, std::vector<double> x,
                std::string colname2, std::vector<double> y);
    DataTable_t(size_t n, size_t m);

    // Modifiers
    void append(std::string colname, std::vector<double> x);
    void resize(Dim_t dim, size_t sz);
    // Acess
    std::vector<double>& at(size_t i);
    double               at(size_t i, size_t j);
    double               at(std::string colname, size_t j);
    // Queries
    size_t getColNumber(std::string colname);
    size_t getNumberRows();
    size_t getNumberCols();
    double getVal(std::string col_out, std::string col_in, double x);

    // Data containers
    std::vector<std::string> cols;
    std::vector< std::vector<double> > data;
};

#endif