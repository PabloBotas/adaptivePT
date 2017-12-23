#include "command_line_parser.hpp"

#include <iostream>
#include <string>
#include <sstream>

#include <exception>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/any.hpp>
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
    po::variables_map vm;
    po::options_description desc("Allowed options");
    try {
        desc.add_options()
        ("help", "Produce this help message.")
        ("skip-cbct",   po::bool_switch(&skip_cbct)->
                        default_value(skip_cbct),
                        "Wheter the calculation on the CBCT should be skiped")
        // Common parameters
        ("patient",     po::value<std::string>(&patient)->
                            required(),
                            "Topas MCAUTO_DICOM.txt file with the plan parameters.")
        ("cbct",        po::value<std::string>(&cbct_file)->
                            required(),
                            "CBCT to adapt the plan to.")
        ("vf",          po::value<std::string>(&vf_file)->
                            default_value(vf_file)->
                            required(),
                            "Vector field file from CT to CBCT. B-Spline format is not supported.")
        ("machine",     po::value<std::string>(&machine)->
                            default_value(machine),
                            "Machine to use during the gPMC simulation or beam model calculation.")
        // Output files
        ("outplan",     po::value<std::string>(&out_plan)->
                            required(),
                            "Output directory to write new gPMC inputs to. "
                            "Will be prepended to any output if they don't contain \'/\'")
        ("optdir",      po::value<std::string>(&out_dir)->
                            required(),
                            "Directory to store optimization and plan evaluation files.")
        ("out_vf",      po::value<std::string>(&data_vf_file)->
                            default_value(data_vf_file)->
                            implicit_value(implicit_data_vf_file),
                            "If the probed values should be written to a file.")
        ("out_shifts",  po::value<std::string>(&data_shifts_file)->
                            default_value(data_shifts_file)->
                            implicit_value(implicit_data_shifts_file),
                            "If the pos-energy shifts should be written to a file.")
        ("traces_ct",   po::value<std::string>(&ct_traces_file)->
                            default_value(ct_traces_file)->
                            implicit_value(implicit_ct_traces_file),
                            "If the traces on the CT volume should be scored to a file.")
        ("traces_cbct", po::value<std::string>(&cbct_traces_file)->
                            default_value(cbct_traces_file)->
                            implicit_value(implicit_cbct_traces_file),
                            "If the traces on the CBCT volume should be scored to a file.")
        ("report",      po::value<std::string>(&report_file)->
                            default_value(report_file)->
                            implicit_value(implicit_report_file),
                            "If a report should be generated. Requires vf_file and data_shifts_file,"
                            " and no geometry only mode")
        // Influence calculations
        ("beam_model",  po::bool_switch(&influence_engine_beam_model)->
                            default_value(influence_engine_beam_model),
                            "If the Dij should be calculated with the native beam model.")
        ("gpmc_dij",    po::bool_switch(&influence_engine_gpmc_dij)->
                            default_value(influence_engine_gpmc_dij),
                            "If the Dij should be calculated with gPMC.")
        ("gpmc_dose",   po::bool_switch(&influence_engine_gpmc_dose)->
                            default_value(influence_engine_gpmc_dose),
                            "If the optimization is restricted to operate in cold/hot areas.")
        ("plan_dose",   po::value<std::string>(&dose_plan_file),
                            "Plan dose distribution, usually from gPMC.")
        ("frac_dose",   po::value<std::string>(&dose_frac_file),
                            "Fraction dose distribution after physical adaptation. This is an "
                            "output file calculated with gPMC and then read for cold/hot spot "
                            "corrections.")
        ("plan_dij",    po::value<std::string>(&dose_frac_file),
                            "Dij file corresponding to the initial plan. Can be interpreted as "
                            "an output file in beam_model mode or as an input file.")
        ("frac_dij",    po::value<std::string>(&dose_frac_file),
                            "Dij file corresponding to the daily geometry plan. It is calculated "
                            "in the daily geometry and after applying the vector field to the "
                            "initial plan. Can be an output from the \"beam_model\" or the "
                            "\"gpmc_dij\" mode")
        ("ct_mask",     po::value<std::string>(&ct_mask_file)->required(),
                            "Binary image of a mask in CT to use during beam model calculation, "
                            "to pass to gPMC for Dij or to evaluate the dose on.")
        // Launchers
        ("opt4D",       po::bool_switch(&launch_opt4D)->
                            default_value(launch_opt4D),
                            "If Opt4D should be launched. If a file destination is not provided as "
                            "well, the default will be set.")
        // Adaptation methods
        // Free positions
        ("free",          po::bool_switch(&FREE_POS_FREE_ENERGY)->
                              default_value(FREE_POS_FREE_ENERGY),
                              "If the positions and energies should move freely.")
        ("rigid-e",       po::bool_switch(&FREE_POS_RIGID_ENERGY)->
                              default_value(FREE_POS_RIGID_ENERGY),
                              "If the positions should move freely and energies rigidly.")
        ("rigid-e-beams", po::bool_switch(&FREE_POS_RIGID_BEAMS_ENERGY)->
                              default_value(FREE_POS_RIGID_BEAMS_ENERGY),
                              "If the positions should move freely and energies rigidly with "
                              "independent fields.")
        // Rigid positions
        ("rigid-pos",         po::bool_switch(&RIGID_POS_FREE_ENERGY)->
                                  default_value(RIGID_POS_FREE_ENERGY),
                                  "If the positions should move rigidly and energies freely.")
        ("rigid",             po::bool_switch(&RIGID_POS_RIGID_ENERGY)->
                                  default_value(RIGID_POS_RIGID_ENERGY),
                                  "If the positions and energies should move rigidly.")
        ("rigid-pos-e-beams", po::bool_switch(&RIGID_POS_RIGID_BEAMS_ENERGY)->
                                  default_value(RIGID_POS_RIGID_BEAMS_ENERGY),
                                  "If the positions should move rigidly and energies rigidly with "
                                  "independent fields.")
        // Rigid per field positions
        ("rigid-pos-beams",   po::bool_switch(&RIGID_BEAMS_POS_FREE_ENERGY)->
                                  default_value(RIGID_BEAMS_POS_FREE_ENERGY),
                                  "If the positions should move rigidly with independent fields and"
                                  " energies freely.")
        ("rigid-pos-beams-e", po::bool_switch(&RIGID_BEAMS_POS_RIGID_ENERGY)->
                                  default_value(RIGID_BEAMS_POS_RIGID_ENERGY),
                                  "If the positions should move rigidly with independent fields and"
                                  " energies rigidly.")
        ("rigid-beams",       po::bool_switch(&RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY)->
                                  default_value(RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY),
                                  "If the positions and energies should move rigidly with "
                                  "independent fields.");

        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (argc == 1)
            throw std::invalid_argument("No input arguments.");
        if (vm.count("help") || argc == 1) {
            std::cerr << desc << std::endl;
            exit(EXIT_SUCCESS);
        }
        po::notify(vm);

        // Check required options
        if (!skip_cbct) {
            if (cbct_file.empty())
                throw std::invalid_argument("the option \'--cbct\' is required but missing");
            if (vf_file.empty())
                throw std::invalid_argument("the option \'--vf\' is required but missing");
        }

        // INFLUENCE OPTIONS
        if(influence_engine_gpmc_dij) {
            influence_opts = Influence_engines_t::GPMC_DIJ;
            if (dose_plan_file.empty())
                throw std::invalid_argument("the plan dose is required when using gPMC for "
                                            "Dij reoptimization");
        }
        else if(influence_engine_gpmc_dose) {
            influence_opts = Influence_engines_t::GPMC_DOSE;
            if (dose_plan_file.empty())
                throw std::invalid_argument("the plan dose is required when using gPMC for "
                                            "gPMC for cold/hot spots correction");
        }
        else
            influence_opts = Influence_engines_t::BEAM_MODEL;
        
        // WARPING OPTIONS
        if (FREE_POS_RIGID_ENERGY)
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
        if ( !report_file.empty() && ( data_shifts_file.empty() || data_vf_file.empty()))
            throw std::invalid_argument("report option needs energies processing and shifts "
                                        "and vf output");

        // CORRECT PATHS
        if (!out_dir.empty() &&
            !data_shifts_file.empty() && data_shifts_file.find('/') == std::string::npos)
            data_shifts_file = out_dir + '/' + data_shifts_file;
        if (!out_dir.empty() &&
            !data_vf_file.empty() && data_vf_file.find('/') == std::string::npos)
            data_vf_file = out_dir + '/' + data_vf_file;
        if (!out_dir.empty() &&
            !ct_traces_file.empty() && ct_traces_file.find('/') == std::string::npos)
            ct_traces_file = out_dir + '/' + ct_traces_file;
        if (!out_dir.empty() &&
            !cbct_traces_file.empty() && cbct_traces_file.find('/') == std::string::npos)
            cbct_traces_file = out_dir + '/' + cbct_traces_file;
    }
    catch(std::exception& e) {
        std::cerr << "ERROR! " << e.what() << std::endl;
        auto_print_map(vm);
        std::cerr << desc << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Parser::print_parameters()
{
    std::cout << "Parsed parameters:\n";
    std::cout << "    - Patient:          " << patient << '\n';
    std::cout << "    - Cone-Beam CT:     " << cbct_file << '\n';
    std::cout << "    - Vector field:     " << vf_file << '\n';
    if (!out_plan.empty())
        std::cout << "    - Out dir:          " << out_plan << '\n';
    if (!data_vf_file.empty())
        std::cout << "    - Out vf:           " << data_vf_file << '\n';
    if (!data_shifts_file.empty())
        std::cout << "    - Out shifts:       " << data_shifts_file << '\n';
    if (!ct_traces_file.empty())
        std::cout << "    - Out CT traces:    " << ct_traces_file << '\n';
    if (!cbct_traces_file.empty())
        std::cout << "    - Out CBCT traces:  " << cbct_traces_file << '\n';
    if (!report_file.empty())
        std::cout << "    - Report:           " << report_file << '\n';
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


void Parser::auto_print_map (po::variables_map vm) {
    std::cout << "The program's inputs were:" << std::endl;
    for (po::variables_map::iterator it = vm.begin(); it != vm.end(); it++) {
        std::cout << "> " << it->first;
        if (((boost::any)it->second.value()).empty()) {
            std::cout << " (empty)";
        }
        if (vm[it->first].defaulted() || it->second.defaulted()) {
            std::cout << " (default)";
        }
        std::cout << " = ";

        bool is_char;
        try {
            boost::any_cast<const char *>(it->second.value());
            is_char = true;
        } catch (const boost::bad_any_cast &) {
            is_char = false;
        }
        bool is_str;
        try {
            boost::any_cast<std::string>(it->second.value());
            is_str = true;
        } catch (const boost::bad_any_cast &) {
            is_str = false;
        }

        if (((boost::any)it->second.value()).type() == typeid(int)) {
            std::cout << vm[it->first].as<int>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(bool)) {
            std::cout << vm[it->first].as<bool>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(float)) {
            std::cout << vm[it->first].as<float>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(double)) {
            std::cout << vm[it->first].as<double>() << std::endl;
        } else if (is_char) {
            std::cout << vm[it->first].as<const char *>() << std::endl;
        } else if (is_str) {
            std::string temp = vm[it->first].as<std::string>();
            if (temp.size()) {
                std::cout << temp << std::endl;
            } else {
                std::cout << "true" << std::endl;
            }
        } else { // Assumes that the only remainder is vector<string>
            try {
                std::vector<std::string> vect = vm[it->first].as<std::vector<std::string> >();
                uint i = 0;
                for (std::vector<std::string>::iterator oit=vect.begin();
                     oit != vect.end(); oit++, ++i) {
                    std::cout << "\r> " << it->first << "[" << i << "]=" << (*oit) << std::endl;
                }
            } catch (const boost::bad_any_cast &) {
                std::cout << "UnknownType(" << ((boost::any)it->second.value()).type().name() << ")" << std::endl;
            }
        }
    }
    std::cout << std::endl;
}
