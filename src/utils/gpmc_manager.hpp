#ifndef __GPMC_MANAGER_HPP__
#define __GPMC_MANAGER_HPP__

#include "enviroment.hpp"
#include "patient_parameters.hpp"
#include "special_types.hpp"
#include "vector4.hpp"

#include <string>
#include <vector>

class Gpmc_manager
{
public:
    Gpmc_manager(Patient_Parameters_t pat_, std::string new_patient,
                 std::string dosefile_,
                 std::string scorer_, std::string outdir_,
                 std::string tramp_dir_,
                 std::vector<std::string> tramp_files_,
                 Vector3_t<float> iso_shift_);
    ~Gpmc_manager();
    void calculate_dij(float min, float max, bool toctgrid, std::vector<std::string> maskfile);
    void calculate_dose(float sf);
    void write_templates(bool add_mask, bool sum_beams);
    void write_dose_files(float sf);
    void write_dij_files(float min_p_spot, float max_p_spot,
                         bool toctgrid, std::vector<std::string> maskfiles);
    void rescale_tramp_weigths(float min_p_spot, float max_p_spot);

    void launch_dij();
    void launch_dose();
    void launch();
    std::string get_total_dose_file();
    std::vector<std::string> get_field_dose_files();
    std::vector<std::string> get_field_dij_files();
    float get_to_Gy_factor();
private:
    Gpmc_manager();
    
    std::string input_file;
    std::string launcher_file;
    std::string patient;
    std::string out_directory;
    std::string machine;
    std::string ct_volume;
    std::vector<RangeShifter_Dims_t> range_shifters;
    std::vector<std::string> tramp_files;
    std::vector<std::string> beam_names;
    Vector3_t<float> iso_shift;

    std::string result_name_root;
    std::string scorer;
    float spot_factor;
    float min_p_per_spot_dij;
    float max_p_per_spot_dij;
    bool to_ct_grid;
    std::vector<std::string> masks;
    std::string trampdir;

    uint nfractions = 30;

    // default file names
    std::string template_launcher_file = std::string(INSTALLATION_PATH) +
                                         "/src/extra/gpmc_launcher_template.sh";
    std::string template_input_file = std::string(INSTALLATION_PATH) +
                                      "/src/extra/gpmc_input_template.in";

};

#endif