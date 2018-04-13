#include "gpmc_manager.hpp"
#include "tramp.hpp"
#include "utils.hpp"

#include <map>

Gpmc_manager::~Gpmc_manager(){}

Gpmc_manager::Gpmc_manager(Patient_Parameters_t plan, std::string patient_,
                           std::string dosefile_, std::string scorer_, std::string outdir_,
                           std::string tramp_dir_, std::vector<std::string> tramp_files_,
                           Vector3_t<float> iso_shift_)
                           : spot_factor(1.f), min_p_per_spot_dij(1e5), max_p_per_spot_dij(1e7)
                             
{
    patient = utils::get_full_path(patient_);
    out_directory = utils::get_full_path(outdir_);
    trampdir = tramp_dir_;
    result_name_root = utils::remove_file_extension(dosefile_);
    machine = plan.machine;
    ct_volume = "./input/ctbinary/ctvolume.dat";
    range_shifters = plan.range_shifters;
    iso_shift = iso_shift_;
    scorer = scorer_;
    beam_names = plan.beam_names;

    for (auto& f: tramp_files_) {
        if (!utils::starts_with_string(f, "/"))
            f = trampdir + "/" + f;
    }
    tramp_files = tramp_files_;

    if (scorer.compare("dosedij") == 0) {
        input_file = out_directory + "/gpmc_dij.in";
        launcher_file = out_directory + "/run_gpmc_dij.sh";
        supress_std = false;
    } else if (scorer.compare("dose") == 0) {
        input_file = out_directory + "/gpmc_dose.in";
        launcher_file = out_directory + "/run_gpmc_dose.sh";
        supress_std = true;
    } else {
        std::cerr << "Unknown scorer passed to gPMC. This should be a compilation error ";
        std::cerr << ", but life is hard sometimes" << std::endl;
    }

}

void Gpmc_manager::calculate_dij(float min_p_spot, float max_p_spot,
                                 bool toctgrid, std::vector<std::string> maskfiles)
{
    write_dij_files(min_p_spot, max_p_spot, toctgrid, maskfiles);
    launch_dij();
}

void Gpmc_manager::write_dij_files(float min_p_spot, float max_p_spot,
                                   bool toctgrid, std::vector<std::string> maskfiles)
{
    to_ct_grid = toctgrid;
    masks = maskfiles;
    min_p_per_spot_dij = min_p_spot;
    max_p_per_spot_dij = max_p_spot;

    rescale_tramp_weigths(min_p_per_spot_dij, max_p_per_spot_dij);
    write_templates(masks.size() > 0, false);
}

void Gpmc_manager::launch_dij()
{
    launch();
}

void Gpmc_manager::calculate_dose(float spot_factor_)
{
    write_dose_files(spot_factor_);
    launch_dose();
}

void Gpmc_manager::write_dose_files(float spot_factor_)
{
    to_ct_grid = true;
    spot_factor = spot_factor_;
    write_templates(false, true);
}

void Gpmc_manager::launch_dose()
{
    launch();
}

void Gpmc_manager::rescale_tramp_weigths(float min_p_spot, float max_p_spot)
{
    // If min_p_spot >= max_p_spot, then it will be set as spot factor
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
    } else {
        spot_factor = min_p_spot;
    }
}


void Gpmc_manager::launch()
{
    std::cout << "Launching gPMC for " << scorer << " scoring!!!" << std::endl;
    std::string cmd("bash " + launcher_file + " 2>&1");
    int return_code{};
    std::string temp = utils::run_command(cmd, &return_code);
    std::cout << temp << std::endl;

    if (return_code != 0) {
        std::cerr << "ERROR running command:" << std::endl;
        std::cerr << cmd << std::endl;
        std::cerr << "Exit code != 0 detected: " << return_code << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Gpmc_manager::write_templates(bool add_mask, bool sum_beams)
{
    std::map<std::string, std::string> replace_map {
        {"PATIENT_DIRECTORY", patient},
        {"SCORER", scorer},
        {"SUPRESS_STD", std::to_string(supress_std)},
        {"SPOT_FACTOR", std::to_string(spot_factor)},
        {"MACHINE_TO_USE", machine},
        {"CT_VOLUME_FILE", ct_volume},
        {"IF_TO_CT_GRID", std::to_string(to_ct_grid)},
        {"OUTPUT_DIR", out_directory},
        {"RESULT_NAME", result_name_root},
        {"TRAMP_DIR", trampdir},
    };
    utils::copy_replace_in_file(template_input_file, input_file, replace_map);

    // Add new range shifters if modified
    // RangeShifter beam1, density1, downstream_face1, thickness1, beam2, density2, downstream_face2, thickness2
    std::string txt("\nRangeShifter");
    bool rs_adapted = false;
    for (size_t i=0; i < range_shifters.size(); ++i) {
        if (range_shifters.at(i).adapted) {
            rs_adapted = true;
            txt += " " + std::to_string(i);
            txt += ", " + std::to_string(range_shifters.at(i).density);
            txt += ", " + std::to_string(range_shifters.at(i).zdown);
            txt += ", " + std::to_string(range_shifters.at(i).thick);
        }
    }
    if (rs_adapted) {
        utils::append_to_file(input_file, txt);
    }

    // Add isocenter shift information
    txt = "\nIsocenterShiftInternalCoords " +
          std::to_string(iso_shift.x) + " " + 
          std::to_string(iso_shift.y) + " " + 
          std::to_string(iso_shift.z);
    if (iso_shift.x != 0.f || iso_shift.y != 0.f || iso_shift.z != 0.f) {
        utils::append_to_file(input_file, txt);
    }

    // Add scoring mask
    if (add_mask) {
        txt = "\nMask ";
        for (size_t i=0; i<masks.size(); ++i) {
            txt += " " + masks.at(i);
            if (i != masks.size()-1)
                txt += ",";
        }
        txt += "\nScopingMask 1";
        utils::append_to_file(input_file, txt);
    }

    replace_map = std::map<std::string, std::string> {
        {"FILE", input_file},
    };
    utils::copy_replace_in_file(template_launcher_file, launcher_file, replace_map);

    if (sum_beams) {
        txt = "\n";
        txt += "fractions=" + std::to_string(nfractions) + "\n";
        txt += "outdir=" + out_directory + "\n";
        txt += "spotfactor=" + std::to_string(spot_factor) + "\n";
        txt += "/opt/gpmc-tools/1.0/sumDoses ${outdir}/total.dose ${outdir}/DoseFrac*.dose\n";
        if (!supress_std) {
            txt += "/opt/gpmc-tools/1.0/sumDosesStd ${outdir}/total.dose_std ${outdir}/DoseFrac*_std\n";
        }
        txt += "nx=$(grep \"nVoxelsX\" ${outdir}/geometry.dat | awk '{print $2}')\n";
        txt += "ny=$(grep \"nVoxelsY\" ${outdir}/geometry.dat | awk '{print $2}')\n";
        txt += "nz=$(grep \"nVoxelsZ\" ${outdir}/geometry.dat | awk '{print $2}')\n";
        txt += "for f in $(ls ${outdir}/total.dose); do\n";
        txt += "    name=$(basename $f)\n";
        txt += "    /opt/gpmc-tools/1.0/gpmc2xio ";
        txt += "${outdir}/ast_${name} ${f} ${nx} ${ny} ${nz} ${spotfactor} ${fractions}\n";
        txt += "done\n";

        utils::append_to_file(launcher_file, txt);
    }

    utils::append_to_file(launcher_file, "\nexit 0\n");
}

std::vector<std::string> Gpmc_manager::get_field_dose_files()
{
    std::vector<std::string> out = beam_names;
    for (auto& str: out) {
        str = out_directory+"/"+result_name_root+"."+str+".dose";
    }
    return out;
}

std::vector<std::string> Gpmc_manager::get_field_dij_files()
{
    std::vector<std::string> out = beam_names;
    for (auto& str: out) {
        str = out_directory+"/"+result_name_root+"."+str+".dij";
    }
    return out;
}

std::string Gpmc_manager::get_total_dose_file()
{
    return out_directory+"/total.dose";
}

float Gpmc_manager::get_to_Gy_factor()
{
    float MeVg2Gy  = 1.6022e-10;
    float RBE       = 1.1;
    float protonsGp = 1.0e9;
    float protonsGpsim = protonsGp / spot_factor;
    return nfractions * RBE * protonsGpsim * MeVg2Gy;
}
