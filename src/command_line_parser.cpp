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
                     "If the traces on the CBCT volume should be scored to a file.");
        // ("layer",   po::value<bool>(&if_per_layer)->default_value(false),
        //             "If the energy should be adapted layer by layer instead of with individual spots")
        // ;

        po::variables_map vm;        
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            exit(EXIT_SUCCESS);
        }

        po::notify(vm);
    }
    catch(std::exception& e) {
        std::cerr << "ERROR! " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!out_dir.empty() && !output_shifts.empty() && output_shifts.find('/') == std::string::npos)
    {
        output_shifts = out_dir + '/' + output_shifts;
    }
    if (!out_dir.empty() && !output_vf.empty() && output_vf.find('/') == std::string::npos)
    {
        output_vf = out_dir + '/' + output_vf;
    }
    if (!out_dir.empty() && !output_ct_traces.empty() && output_ct_traces.find('/') == std::string::npos)
    {
        output_ct_traces = out_dir + '/' + output_ct_traces;
    }
    if (!out_dir.empty() && !output_cbct_traces.empty() && output_cbct_traces.find('/') == std::string::npos)
    {
        output_cbct_traces = out_dir + '/' + output_cbct_traces;
    }
}

// bool Parser::check_both_or_none(po::variables_map vm, std::string arg1, std::string arg2)
// {
//     if(vm.count(arg1) != vm.count(arg2))
//     {
//         std::stringstream errMsg;
//         errMsg << arg1 << " and " << arg2 << " must BOTH be specified or omitted." << std::endl;
//         throw std::invalid_argument( errMsg.str() );
//     }
//     bool both = vm.count(arg1);
//     return both;
// }

// void Parser::check_vector_sizes(po::variables_map vm, std::string arg1, std::string arg2)
// {
//     std::vector<std::string> vec1 = vm[arg1].as< std::vector<std::string> >();
//     std::vector<std::string> vec2 = vm[arg2].as< std::vector<std::string> >();
//     if(vec1.size() != vec2.size())
//     {
//         for (const auto i: vec1)
//             std::cout << i << ' ';
//         std::cout << std::endl;

//         for (const auto i: vec2)
//             std::cout << i << ' ';
//         std::cout << std::endl;

//         std::stringstream errMsg;
//         errMsg << arg1 << " and " << arg2 << " must have the same number of elements." << std::endl;
//         throw std::invalid_argument( errMsg.str() );
//     }
// }


void Parser::print_parameters()
{
    std::cout << "Parsed parameters:\n";
    std::cout << "    - Patient:      " << patient << '\n';
    std::cout << "    - Cone-Beam CT: " << cbct_file << '\n';
    std::cout << "    - Vector field: " << vf_file << '\n';
    if (!out_dir.empty())
        std::cout << "    - Out dir:      " << out_dir << '\n';
}

