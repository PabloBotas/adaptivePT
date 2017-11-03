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

    // General control
    bool skip_cbct;
    // Common parameters
    std::string patient;
    std::string cbct_file;
    std::string vf_file;
    // Output files
    std::string out_dir;
    std::string output_shifts;
    std::string output_vf;
    std::string output_ct_traces;
    std::string report;
    std::string output_cbct_traces;
    std::string output_opt4D_files;
    bool launch_opt4D;
    std::string ct_target_file;
    // Adaptation methods
    Warp_opts_t warp_opts;

    void print_parameters();

private:
    bool FREE_POS_FREE_ENERGY;
    bool FREE_POS_RIGID_ENERGY;
    bool FREE_POS_RIGID_BEAMS_ENERGY;

    bool RIGID_POS_FREE_ENERGY;
    bool RIGID_POS_RIGID_ENERGY;
    bool RIGID_POS_RIGID_BEAMS_ENERGY;

    bool RIGID_BEAMS_POS_FREE_ENERGY;
    bool RIGID_BEAMS_POS_RIGID_ENERGY;
    bool RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY;

    // Default and implicit values (implicit values act as default when a non-
    // required parameter is passed with no arguments)
    bool default_skip_cbct                  = false;
    std::string default_cbct_file           = std::string();
    std::string default_vf_file             = std::string();
    std::string default_output_vf           = std::string();
    std::string implicit_output_vf          = std::string("vf.dat");
    std::string default_output_shifts       = std::string();
    std::string implicit_output_shifts      = std::string("shifts.dat");
    std::string default_output_ct_traces    = std::string();
    std::string implicit_output_ct_traces   = std::string("ct_traces.mhd");
    std::string default_output_cbct_traces  = std::string();
    std::string implicit_output_cbct_traces = std::string("cbct_traces.pdf");
    std::string default_report              = std::string();
    std::string implicit_report             = std::string("adapt_report.pdf");
    std::string default_output_opt4D_files  = std::string();
    std::string implicit_output_opt4D_files = std::string("opt4D_files");
    bool default_launch_opt4D                       = false;
    bool default_FREE_POS_FREE_ENERGY               = false;
    bool default_FREE_POS_RIGID_ENERGY              = false;
    bool default_FREE_POS_RIGID_BEAMS_ENERGY        = false;
    bool default_RIGID_POS_FREE_ENERGY              = false;
    bool default_RIGID_POS_RIGID_ENERGY             = false;
    bool default_RIGID_POS_RIGID_BEAMS_ENERGY       = false;
    bool default_RIGID_BEAMS_POS_FREE_ENERGY        = false;
    bool default_RIGID_BEAMS_POS_RIGID_ENERGY       = false;
    bool default_RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY = false;


    void process_command_line(int argc, char** argv);
    // bool check_both_or_none(po::variables_map vm, std::string arg1, std::string arg2);
    // void check_vector_sizes(po::variables_map vm, std::string arg1, std::string arg2);
};

#endif
