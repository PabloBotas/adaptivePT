#ifndef __CMD_LINE_PARSER__
#define __CMD_LINE_PARSER__

#include <string>
#include <vector>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

class Parser
{
public:
    Parser(int argc, char** argv);
    ~Parser();

    std::string patient;
    std::string cbct_file;
    std::string vf_file;
    std::string out_dir;
    std::string output_shifts;
    std::string output_vf;
    std::string output_ct_traces;
    bool ct_traces_individual;
    std::string output_cbct_traces;
    bool if_per_layer;
    bool no_energy;

    void print_parameters();

private:
    void process_command_line(int argc, char** argv);
    // bool check_both_or_none(po::variables_map vm, std::string arg1, std::string arg2);
    // void check_vector_sizes(po::variables_map vm, std::string arg1, std::string arg2);
};

#endif
