#ifndef __CMD_LINE_PARSER__
#define __CMD_LINE_PARSER__

#include "program_options.hpp"
#include <string>
#include <vector>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

class Parser
{
public:
    Parser(int argc, char** argv);
    ~Parser();

    // General control and debugging tools
    bool skip_cbct = false;
    // Common parameters
    std::string patient;
    std::string cbct_file;
    std::string vf_file;
    std::string machine = "topasmediumspots";
    // Output files
    std::string out_plan;
    std::string out_dir;
    std::string data_shifts_file;
    std::string data_vf_file;
    std::string ct_traces_file;
    std::string cbct_traces_file;
    std::string report_file;
    // Launchers
    bool launch_opt4D = false;
    bool launch_adapt_simulation = false;
    // std::string ct_mask_file;
    std::vector<std::string> scoring_mask_files;
    std::vector<std::string> v_field_mask_files;
    // Adaptation method
    Adapt_methods_t adapt_method;
    Adapt_constraints_t adapt_constraints;
    RShifter_steps_t rshifter_steps;
    // Adaptation data files
    std::string dose_plan_file;
    std::string dose_frac_file;
    std::string dij_plan_file;
    std::string dij_frac_file;

    void print_inputs();

private:
    std::string adapt_method_str;
    std::string adapt_method_geometric_str  = "geometric";
    std::string adapt_method_beam_model_str = "beam_model";
    std::string adapt_method_gpmc_dij_str   = "gpmc_dij";
    std::string adapt_method_cold_spots_str = "cold_spots";
    std::vector<std::string> constraint_vec;
    std::string adapt_constraint_free_str      = "free";
    std::string adapt_constraint_v_rs_str      = "virt_range_shifter";
    std::string adapt_constraint_rs_str        = "range_shifter";
    std::string adapt_constraint_iso_str       = "iso_shift";
    std::string adapt_constraint_iso_field_str = "iso_shift_field";
    std::vector<std::string> rshifter_steps_vec;
    std::string rshifter_steps_free    = "free";
    std::string rshifter_steps_half_cm = "half_cm";
    std::string rshifter_steps_cm      = "cm";
    std::string rshifter_steps_mgh     = "mgh";

    // std::string scoring_mask_files_str;
    // std::string v_field_mask_files_str;

    // SHIFTS CONSTRAINS VARIABLES
    bool free = false;
    bool v_range_shifter = false;
    bool range_shifter = false;
    bool iso_shift = false;
    bool iso_shift_field = false;

    // Influence engines variables
    bool influence_engine_beam_model = false;
    bool influence_engine_gpmc_dij   = false;
    bool influence_engine_gpmc_dose  = false;

    // Default and implicit values (implicit values act as default when a non-
    // required parameter is passed with no arguments)
    std::string implicit_data_vf_file     = std::string("vf.dat");
    std::string implicit_data_shifts_file = std::string("shifts.dat");
    std::string implicit_ct_traces_file   = std::string("ct_traces.mhd");
    std::string implicit_cbct_traces_file = std::string("cbct_traces.mhd");
    std::string implicit_report_file      = std::string("adapt_report.pdf");
    std::string default_dose_frac_file    = std::string("DoseFrac");
    std::string default_dij_frac_file     = std::string("DijFrac");
    

    void process_command_line(int argc, char** argv);
    std::vector<std::string> get_string_items(const std::string& str, char sep=',') {
        std::vector<std::string> out;
        if (!str.empty()) {
            std::stringstream stream(str);
            while (stream.good()) {
                std::string temp;
                getline(stream, temp, sep);
                temp = std::regex_replace(temp, std::regex("^ +| +$"), "");
                out.push_back(temp);
            }
        }
        return out;
    };
    // bool check_both_or_none(po::variables_map vm, std::string arg1, std::string arg2);
    // void check_vector_sizes(po::variables_map vm, std::string arg1, std::string arg2);
    void query_terminal_width();
    int terminal_width;

    po::variables_map vm;
};

#endif
