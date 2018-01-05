#ifndef __OPT4D_MANAGER_HPP__
#define __OPT4D_MANAGER_HPP__

#include "enviroment.hpp"
#include "patient_parameters.hpp"
#include "vector4.hpp"

#include <string>
#include <vector>

class Gpmc_manager
{
public:
    Gpmc_manager(Patient_Parameters_t pat_, std::string dosefile_,
                 std::string scorer_, std::string trampdir_,
                 std::vector<std::string> tramp_files_);
    ~Gpmc_manager();
    void write_templates();
    void rescale_tramp_weigths(float min_p_spot, float max_p_spot);
    void calculate_dij(float min, float max, bool toctgrid, std::string maskfile);
    void launch();
private:
    Gpmc_manager();
    
    std::string input_file;
    std::string launcher_file;
    std::string pat_directory;
    std::string out_directory;
    std::string machine;
    std::string ct_volume;
    std::vector<std::string> tramp_files;

    std::string result_name;
    std::string dose_file;
    std::string scorer;
    float spot_factor;
    float min_p_per_spot_dij;
    float max_p_per_spot_dij;
    bool to_ct_grid;
    std::string mask;
    std::string trampdir;
    
    // default file names
    std::string template_launcher_file = std::string(INSTALLATION_PATH) +
                                         "/src/extra/gpmc_launcher_template.sh";
    std::string template_input_file = std::string(INSTALLATION_PATH) +
                                      "/src/extra/gpmc_input_template.in";

};

#endif