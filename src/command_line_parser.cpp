#include "command_line_parser.hpp"

#include <iostream>
#include <string>
#include <sstream>

#include <exception>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

Parser::Parser(int argc, char** argv)
{
    process_command_line(argc, argv);
}

Parser::~Parser()
{
}

void Parser::process_command_line(int argc, char** argv)
{
    try 
    {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "Produce this help message.")
        // Common parameters
        ("patient", po::value<std::string>(&patient)->required(),
                    "Topas MCAUTO_DICOM.txt file with the plan parameters.")
        ("cbct",    po::value<std::string>(&cbct_file)->required(),
                    "CBCT to adapt the plan to.")
        ("vf",      po::value<std::string>(&vf_file)->required(),
                    "Vector field file from CT to CBCT. B-Spline format is not supported.")
        ("outdir",  po::value<std::string>(&out_dir)->required(),
                    "Output directory to write results to. Will be prepended to any output if they don't contain \'/\'")
        ("no_energy", po::bool_switch(&no_energy)->default_value(false),
                    "If only the geometry should be adapted.")
        ("output_vf", po::value<std::string>(&output_vf)->default_value(std::string()),
                    "If the probed values should be written to a file.")
        ("output_shifts", po::value<std::string>(&output_shifts)->default_value(std::string()),
                    "If the pos-energy shifts should be written to a file.")
        ("output_ct_traces", po::value<std::string>(&output_ct_traces)->default_value(std::string()),
                    "If the traces on the CT volume should be scored to a file.")
        ("output_cbct_traces", po::value<std::string>(&output_cbct_traces)->default_value(std::string()),
                    "If the traces on the CBCT volume should be scored to a file.")
        ("ct_traces_individual", po::bool_switch(&ct_traces_individual)->default_value(false),
                    "If the traces should be scored individually.")
        ("report", po::bool_switch(&report)->default_value(false),
                    "If a report should be generated. Requires output_vf and output_shifts");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            std::cerr << desc << std::endl;
            exit(EXIT_SUCCESS);
        }
        po::notify(vm);
    
        if ( ct_traces_individual ) {
            if ( !output_ct_traces.empty() && !output_cbct_traces.empty() ) {
                std::cerr << "ERROR! ct_traces_individual option needs output_ct_traces or output_cbct_traces." << std::endl;
                std::cerr << desc << std::endl;
                exit(EXIT_SUCCESS);
            }
        }
        if ( report ) {
            if ( !output_shifts.empty() || !output_vf.empty() || no_energy) {
                std::cerr << "ERROR! report option needs energies processing and shifts and vf output." << std::endl;
                std::cerr << desc << std::endl;
                exit(EXIT_SUCCESS);
            }
        }

        if (!out_dir.empty() && !output_shifts.empty() && output_shifts.find('/') == std::string::npos)
            output_shifts = out_dir + '/' + output_shifts;
        if (!out_dir.empty() && !output_vf.empty() && output_vf.find('/') == std::string::npos)
            output_vf = out_dir + '/' + output_vf;
        if (!out_dir.empty() && !output_ct_traces.empty() && output_ct_traces.find('/') == std::string::npos)
            output_ct_traces = out_dir + '/' + output_ct_traces;
        if (!out_dir.empty() && !output_cbct_traces.empty() && output_cbct_traces.find('/') == std::string::npos)
            output_cbct_traces = out_dir + '/' + output_cbct_traces;
    }
    catch(std::exception& e) {
        std::cerr << "ERROR! " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Parser::print_parameters()
{
    std::cout << "Parsed parameters:\n";
    std::cout << "    - Patient:          " << patient << '\n';
    std::cout << "    - Cone-Beam CT:     " << cbct_file << '\n';
    std::cout << "    - Vector field:     " << vf_file << '\n';
    if (!out_dir.empty())
        std::cout << "    - Out dir:          " << out_dir << '\n';
    if (no_energy)
        std::cout << "    - No Energy" << '\n';
    if (!output_vf.empty())
        std::cout << "    - Out vf:           " << output_vf << '\n';
    if (!output_shifts.empty())
        std::cout << "    - Out shifts:       " << output_shifts << '\n';
    if (!output_ct_traces.empty())
        std::cout << "    - Out CT traces:    " << output_ct_traces << '\n';
    if (!output_cbct_traces.empty())
        std::cout << "    - Out CBCT traces:  " << output_cbct_traces << '\n';
    if (ct_traces_individual)
        std::cout << "    - Traces Individual" << '\n';
    if (report)
        std::cout << "    - Report" << '\n';
}
