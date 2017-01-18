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

Parser::Parser()
{
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
                    "Topas MCAUTO_DICOM.txt file with the plan parameters")
        ("cbct",    po::value<std::string>(&cbct_file)->required(),
                    "CBCT to adapt the plan to")
        ("layer",   po::value<bool>(&if_per_layer)->default_value(false),
                    "If the energy should be adapted layer by layer instead of with individual spots")
        ;

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
    // std::cout << "    - tramp_in_file: " << tramp_in_file << '\n';
    // std::cout << "    - dvh_file:      " << dvh_file << '\n';
    // std::cout << "    - target_dose:   " << target_dose << " Gy" << '\n';
    // if(perform_Bwf_scaling)
    // {
    //     std::cout << "    - bwf_in_file:   " << bwf_in_file << '\n';
    //     std::cout << "    - bwf_out_file:  " << bwf_out_file << '\n';
    // }
    // if(perform_operate_volume)
    // {
    //     std::cout << "    - dose_in_files (" << dose_in_file.size() << "): \n";
    //     for (auto i : dose_in_file)
    //         std::cout << i << ' ';
    //     std::cout << '\n';
    //     std::cout << "    - dose_out_file (" << dose_out_file.size() << "): \n";
    //     for (auto i : dose_out_file)
    //         std::cout << i << ' ';
    //     std::cout << '\n';
    //     std::cout << "    - data_type:     " << data_type.getTypeName() << std::endl;
    // }
}

