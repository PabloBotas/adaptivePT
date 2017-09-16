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
        ("patient",     po::value<std::string>(&patient)->required(),
                        "Topas MCAUTO_DICOM.txt file with the plan parameters.")
        ("cbct",        po::value<std::string>(&cbct_file)->required(),
                        "CBCT to adapt the plan to.")
        ("vf",          po::value<std::string>(&vf_file)->required(),
                        "Vector field file from CT to CBCT. B-Spline format is not supported.")
        ("dose",        po::value<std::string>(&dose_file),
                        "Dose of plan evaluated in CBCT.")
        ("presc",       po::value<double>(&dose_presc),
                        "Target dose prescription.")
        ("spot_factor", po::value<double>(&spot_factor),
                        "Spot factor used in simulation.")
        // Output files
        ("outdir",      po::value<std::string>(&out_dir)->required(),
                        "Output directory to write results to. Will be prepended to any output if they don't contain \'/\'")
        ("out_vf",      po::value<std::string>(&output_vf)->default_value(std::string())->implicit_value("vf.dat"),
                        "If the probed values should be written to a file.")
        ("out_shifts",  po::value<std::string>(&output_shifts)->default_value(std::string())->implicit_value("shifts.dat"),
                        "If the pos-energy shifts should be written to a file.")
        ("traces_ct",   po::value<std::string>(&output_ct_traces)->default_value(std::string())->implicit_value("ct_traces.mhd"),
                        "If the traces on the CT volume should be scored to a file.")
        ("traces_cbct", po::value<std::string>(&output_cbct_traces)->default_value(std::string())->implicit_value("cbct_traces.pdf"),
                        "If the traces on the CBCT volume should be scored to a file.")
        ("report",      po::value<std::string>(&report)->default_value(std::string())->implicit_value("adapt_report.pdf"),
                        "If a report should be generated. Requires output_vf and output_shifts, and no geometry only mode")
        // Adaptation methods
        // Free positions
        ("free",          po::bool_switch(&FREE_POS_FREE_ENERGY)->default_value(false),
                          "If the positions and energies should move freely.")
        ("rigid-e",       po::bool_switch(&FREE_POS_RIGID_ENERGY)->default_value(false),
                          "If the positions should move freely and energies rigidly.")
        ("rigid-e-beams", po::bool_switch(&FREE_POS_RIGID_BEAMS_ENERGY)->default_value(false),
                          "If the positions should move freely and energies rigidly with independent fields.")
        // Rigid positions
        ("rigid-pos",         po::bool_switch(&RIGID_POS_FREE_ENERGY)->default_value(false),
                              "If the positions should move rigidly and energies freely.")
        ("rigid",             po::bool_switch(&RIGID_POS_RIGID_ENERGY)->default_value(false),
                              "If the positions and energies should move rigidly.")
        ("rigid-pos-e-beams", po::bool_switch(&RIGID_POS_RIGID_BEAMS_ENERGY)->default_value(false),
                              "If the positions should move rigidly and energies rigidly with independent fields.")
        // Rigid per field positions
        ("rigid-pos-beams",   po::bool_switch(&RIGID_BEAMS_POS_FREE_ENERGY)->default_value(false),
                              "If the positions should move rigidly with independent fields and energies freely.")
        ("rigid-pos-beams-e", po::bool_switch(&RIGID_BEAMS_POS_RIGID_ENERGY)->default_value(false),
                              "If the positions should move rigidly with independent fields and energies rigidly.")
        ("rigid-beams",       po::bool_switch(&RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY)->default_value(false),
                              "If the positions and energies should move rigidly with independent fields.");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            std::cerr << desc << std::endl;
            exit(EXIT_SUCCESS);
        }
        po::notify(vm);
    
        // WARPING OPTIONS
        if(FREE_POS_RIGID_ENERGY)
            warp_opts = Warp_opts_t::FREE_POS_RIGID_ENERGY;
        else if(FREE_POS_RIGID_BEAMS_ENERGY)
            warp_opts = Warp_opts_t::FREE_POS_RIGID_BEAMS_ENERGY;
       else if(RIGID_POS_FREE_ENERGY)
            warp_opts = Warp_opts_t::RIGID_POS_FREE_ENERGY;
        else if(RIGID_POS_RIGID_ENERGY)
            warp_opts = Warp_opts_t::RIGID_POS_RIGID_ENERGY;
        else if(RIGID_POS_RIGID_BEAMS_ENERGY)
            warp_opts = Warp_opts_t::RIGID_POS_RIGID_BEAMS_ENERGY;
        else if(RIGID_BEAMS_POS_FREE_ENERGY)
            warp_opts = Warp_opts_t::RIGID_BEAMS_POS_FREE_ENERGY;
        else if(RIGID_BEAMS_POS_RIGID_ENERGY)
            warp_opts = Warp_opts_t::RIGID_BEAMS_POS_RIGID_ENERGY;
        else if(RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY)
            warp_opts = Warp_opts_t::RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY;
        else
            warp_opts = Warp_opts_t::FREE_POS_FREE_ENERGY;


        // REPORTING
        if ( !report.empty() ) {
            if ( output_shifts.empty() || output_vf.empty()) {
                std::cerr << "ERROR! report option needs energies processing and shifts and vf output." << std::endl;
                print_parameters();
                std::cerr << desc << std::endl;
                exit(EXIT_SUCCESS);
            }
        }

        // CORRECT PATHS
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
    if (!output_vf.empty())
        std::cout << "    - Out vf:           " << output_vf << '\n';
    if (!output_shifts.empty())
        std::cout << "    - Out shifts:       " << output_shifts << '\n';
    if (!output_ct_traces.empty())
        std::cout << "    - Out CT traces:    " << output_ct_traces << '\n';
    if (!output_cbct_traces.empty())
        std::cout << "    - Out CBCT traces:  " << output_cbct_traces << '\n';
    if (!report.empty())
        std::cout << "    - Report:           " << report << '\n';
    // Warping options
    if (FREE_POS_FREE_ENERGY)
        std::cout << "    - Warping: Free position and energy" << '\n';
    else if(FREE_POS_RIGID_ENERGY)
        std::cout << "    - Warping: Free position and rigid energy" << '\n';
    else if(FREE_POS_RIGID_BEAMS_ENERGY)
        std::cout << "    - Warping: Free position and rigid fields energy" << '\n';
    else if(RIGID_POS_FREE_ENERGY)
        std::cout << "    - Warping: Rigid position and free energy" << '\n';
    else if(RIGID_POS_RIGID_ENERGY)
        std::cout << "    - Warping: Rigid position and energy" << '\n';
    else if(RIGID_POS_RIGID_BEAMS_ENERGY)
        std::cout << "    - Warping: Rigid position and rigid fields energy" << '\n';
    else if(RIGID_BEAMS_POS_FREE_ENERGY)
        std::cout << "    - Warping: Rigid fields position and free energy" << '\n';
    else if(RIGID_BEAMS_POS_RIGID_ENERGY)
        std::cout << "    - Warping: Rigid fields position and rigid energy" << '\n';
    else if(RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY)
        std::cout << "    - Warping: Rigid fields position and energy" << '\n';
}
