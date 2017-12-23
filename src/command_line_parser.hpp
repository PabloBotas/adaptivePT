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
    bool launch_opt4D = false;
    std::string ct_mask_file;
    // Influence engines
    Influence_engines_t influence_opts;
    std::string dose_plan_file;
    std::string dose_frac_file;
    // std::string gpmc_dij_plan_file;
    std::string gpmc_dij_frac_file;
    std::string beam_model_dij_plan_file;
    std::string beam_model_dij_frac_file;
    // Adaptation methods
    Warp_opts_t warp_opts;

    void print_parameters();
    void auto_print_map(po::variables_map vm);

private:
    // SHIFTS CONSTRAINS VARIABLES
    bool FREE_POS_FREE_ENERGY               = false;
    bool FREE_POS_RIGID_ENERGY              = false;
    bool FREE_POS_RIGID_BEAMS_ENERGY        = false;
    bool RIGID_POS_FREE_ENERGY              = false;
    bool RIGID_POS_RIGID_ENERGY             = false;
    bool RIGID_POS_RIGID_BEAMS_ENERGY       = false;
    bool RIGID_BEAMS_POS_FREE_ENERGY        = false;
    bool RIGID_BEAMS_POS_RIGID_ENERGY       = false;
    bool RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY = false;

    // Influence engines variables
    bool influence_engine_beam_model = false;
    bool influence_engine_gpmc_dij   = false;
    bool influence_engine_gpmc_dose  = false;

    // Default and implicit values (implicit values act as default when a non-
    // required parameter is passed with no arguments)
    std::string implicit_data_vf_file     = std::string("vf.dat");
    std::string implicit_data_shifts_file = std::string("shifts.dat");
    std::string implicit_ct_traces_file   = std::string("ct_traces.mhd");
    std::string implicit_cbct_traces_file = std::string("cbct_traces.pdf");
    std::string implicit_report_file      = std::string("adapt_report.pdf");
    

    void process_command_line(int argc, char** argv);
    // bool check_both_or_none(po::variables_map vm, std::string arg1, std::string arg2);
    // void check_vector_sizes(po::variables_map vm, std::string arg1, std::string arg2);
};

#endif
