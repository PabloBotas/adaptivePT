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
    DataTable_t(std::string colname, std::vector<float> x);
    DataTable_t(std::string colname1, std::vector<float> x,
                std::string colname2, std::vector<float> y);
    DataTable_t(size_t n, size_t m);

    // Modifiers
    void append(std::string colname, std::vector<float> x);
    void resize(Dim_t dim, size_t sz);
    // Acess
    std::vector<float>& at(size_t i);
    float               at(size_t i, size_t j);
    float               at(std::string colname, size_t j);
    // Queries
    size_t getColNumber(std::string colname);
    size_t getNumberRows();
    size_t getNumberCols();
    float getVal(std::string col_out, std::string col_in, float x);

    // Data containers
    std::vector<std::string> cols;
    std::vector< std::vector<float>> data;
};

#endif