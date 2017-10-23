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
    std::string dose_file;
    double dose_presc;
    double spot_factor;
    // Output files
    std::string out_dir;
    std::string output_shifts;
    std::string output_vf;
    std::string output_ct_traces;
    std::string report;
    std::string output_cbct_traces;
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

    void process_command_line(int argc, char** argv);
    // bool check_both_or_none(po::variables_map vm, std::string arg1, std::string arg2);
    // void check_vector_sizes(po::variables_map vm, std::string arg1, std::string arg2);
};

#endif
