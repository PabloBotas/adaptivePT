#ifndef __OPT4D_MANAGER_HPP__
#define __OPT4D_MANAGER_HPP__

#include "enviroment.hpp"
#include "vector4.hpp"

#include <string>
#include <vector>

class Opt4D_manager
{
public:
    Opt4D_manager(std::string outdir);
    ~Opt4D_manager();
    void populate_directory(const Array4<double>& influence1,
                            const Array4<double>& influence2);
    void launch_optimization();
    std::vector<double> get_weight_scaling();
private:
    Opt4D_manager();
    void read_bwf_file();
    void set_write_reference_influence(const Array4<double>& influence);
    void write_templates();
    void write_dif();
    void write_vv();
    void write_dij(const Array4<double>& influence);
    
    unsigned int n;
    std::string out_directory;
    double min_average_constrain;
    std::vector<double> weight_scaling;

    // default file names
    std::string launcher_file;
    std::string launcher_file_base = "opt4D_launcher.sh";
    std::string template_launcher_file = std::string(INSTALLATION_PATH) +
                                         "/src/extra/opt4D_launcher_template.sh";
    std::string plan_file;
    std::string plan_file_base = "opt4D_planfile.pln";
    std::string template_plan_file = std::string(INSTALLATION_PATH) +
                                     "/src/extra/opt4D_planfile_template.pln";
    std::string dif_file;
    std::string dif_file_base = "dimensions.dif";
    std::string vv_file;
    std::string vv_file_base = "structures.vv";
    std::string dij_file;
    std::string dij_file_base = "beam_1.dij";
    std::string reference_file;
    std::string reference_file_base = "reference_influence.dat";
    std::string bwf_file;
    std::string bwf_file_base = "beam_1.bwf";

};

#endif