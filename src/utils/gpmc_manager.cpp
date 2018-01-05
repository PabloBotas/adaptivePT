#include "gpmc_manager.hpp"
#include "tramp.hpp"
#include "utils.hpp"

#include <map>

Gpmc_manager::~Gpmc_manager(){}

Gpmc_manager::Gpmc_manager(Patient_Parameters_t pat_, std::string dosefile_,
                           std::string scorer_, std::string outdir_,
                           std::vector<std::string> tramp_files_)
                           : spot_factor(1.f), min_p_per_spot_dij(1e5), max_p_per_spot_dij(1e7)
                             
{
    pat_directory = utils::get_full_path(pat_.patient_dir);
    out_directory = utils::get_full_path(outdir_);
    trampdir = utils::get_full_path(outdir_);
    dose_file = dosefile_;
    result_name = utils::remove_file_extension(dose_file);
    machine = pat_.machine;
    ct_volume = utils::remove_file_extension(pat_.planning_ct_file);
    scorer = scorer_;

    for (auto& f: tramp_files_) {
        if (!utils::starts_with_string(f, "/"))
            f = trampdir + "/" + f;
    }
    tramp_files = tramp_files_;

    if (scorer.compare("dosedij") == 0) {
        input_file = out_directory + "/gpmc_dij.in";
        launcher_file = out_directory + "/run_gpmc_dij.sh";
    } else if (scorer.compare("dose") == 0) {
        input_file = out_directory + "/gpmc_dose.in";
        launcher_file = out_directory + "/run_gpmc_dose.sh";
    }

}


void Gpmc_manager::calculate_dij(float min_p_spot, float max_p_spot,
                                 bool toctgrid, std::string maskfile)
{
    to_ct_grid = toctgrid;
    mask = maskfile;
    min_p_per_spot_dij = min_p_spot;
    max_p_per_spot_dij = max_p_spot;

    write_templates();
    rescale_tramp_weigths(min_p_per_spot_dij, max_p_per_spot_dij);
    launch();
}


void Gpmc_manager::rescale_tramp_weigths(float min_p_spot, float max_p_spot)
{
    if (min_p_spot < max_p_spot) {
        float min_weight = 10000000;
        float max_weight = 0.f;
        for (auto f: tramp_files) {
            Tramp_t tramp(f);
            std::vector<float> weights = tramp.get_weights();
            for (auto w: weights) {
                min_weight = std::min(min_weight, w);
                max_weight = std::max(max_weight, w);
            }
        }
        for (auto f: tramp_files) {
            Tramp_t tramp(f);
            std::vector<float> weights = tramp.get_weights();
            for (auto& w: weights) {
                w = min_p_spot + (max_p_spot-min_p_spot) / (max_weight-min_weight) * (w-min_weight);
            }
            tramp.set_weights(weights);
            tramp.to_file(f);
        }
    }
}


void Gpmc_manager::launch()
{
    std::cerr << "LAUNCHING GPMC!!! not implemented!" << std::endl;
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

    replace_map = std::map<std::string, std::string> {
        {"FILE", input_file},
    };
    utils::copy_file(template_launcher_file, launcher_file, replace_map);
}
