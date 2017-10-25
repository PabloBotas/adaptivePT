#ifndef __OPT4D_MANAGER_HPP__
#define __OPT4D_MANAGER_HPP__

#include "vector4.hpp"

#include <string>

class Opt4D_manager
{
public:
    Opt4D_manager(std::string outdir);
    ~Opt4D_manager();
    void populate_directory(const Array4<double>& influence_ct,
                            const Array4<double>& influence_cbct);
    void launch_optimization();
private:
    Opt4D_manager();
    void write_templates();
    void write_dif();
    void write_vv();
    void write_dij(const Array4<double>& influence_ct);
    void write_reference_influence(const Array4<double>& influence_cbct);
    std::string out_directory;
    unsigned int n;
};

#endif