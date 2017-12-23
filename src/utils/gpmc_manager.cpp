#include "gpmc_manager.hpp"
#include "utils.hpp"

#include <map>

Gpmc_manager::~Gpmc_manager(){}

Gpmc_manager::Gpmc_manager(Patient_Parameters_t pat_, std::string dosefile_,
                           std::string scorer_, std::string trampdir_)
{
    pat_directory = pat_.patient_dir;
    trampdir = trampdir_;
    machine = pat_.machine;
    ct_volume = pat_.planning_ct_file;
    scorer = scorer_;

    std::string temp = utils::get_parent_directory(dosefile_);
    if (utils::starts_with_string(temp, "/"))
        out_directory = temp;
    else
        out_directory = pat_directory + "/" + temp;

    if (scorer.compare("dosedij") == 0) {
        input_file = "gpmc_dij.in";
        launcher_file = "run_gpmc_dij.sh";
    } else if (scorer.compare("dose") == 0) {
        input_file = "gpmc_dose.in";
        launcher_file = "run_gpmc_dose.in";
    }

    dose_file = dosefile_;
    result_name = utils::remove_file_extension(dose_file);
}


void Gpmc_manager::calculate_dij(float spotfactor, bool toctgrid, std::string maskfile)
{
    spot_factor = spotfactor;
    to_ct_grid = toctgrid;
    mask = maskfile;

    write_templates();
    launch();
}


void Gpmc_manager::launch()
{
    std::cout << "LAUNCHING GPMC!!! not implemented!" << std::endl;
    exit(EXIT_SUCCESS);
}


void Gpmc_manager::write_templates()
{
    std::map<std::string, std::string> replace_map {
        {"PATIENT_DIRECTORY", pat_directory},
        {"SCORER", scorer},
        {"SPOT_FACTOR", std::to_string(spot_factor)},
        {"MACHINE_TO_USE", machine},
        {"CT_VOLUME_FILE", ct_volume},
        {"IF_TO_CT_GRID", std::to_string(to_ct_grid)},
        {"OUTPUT_DIR", out_directory},
        {"RESULT_NAME", result_name},
        {"TRAMP_DIR", trampdir},
    };

    utils::copy_file(template_input_file, input_file, replace_map);

    replace_map = std::map<std::string, std::string>{
        {"FILE", pat_directory},
    };
    utils::copy_file(template_launcher_file, launcher_file, replace_map);
}
