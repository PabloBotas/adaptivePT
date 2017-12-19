#ifndef __OPT4D_MANAGER_HPP__
#define __OPT4D_MANAGER_HPP__

#include "enviroment.hpp"
#include "vector4.hpp"

#include <string>
#include <vector>

class Gpmc_manager
{
public:
    Gpmc_manager(std::string outdir);
    ~Gpmc_manager();
    void populate_directory(const Array4<float>& influence1,
                            const Array4<float>& influence2);
    void launch_calculation();
private:
    Gpmc_manager();
    
    std::string out_directory;

    // default file names
    std::string launcher_file;
    std::string launcher_file_base = "opt4D_launcher.sh";
    std::string template_launcher_file = std::string(INSTALLATION_PATH) +
                                         "/src/extra/gpmc_launcher_template.sh";
    std::string plan_file;
    std::string plan_file_base = "opt4D_planfile.pln";
    std::string template_plan_file = std::string(INSTALLATION_PATH) +
                                     "/src/extra/gpmc_input_template.in";
};

#endif