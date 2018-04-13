#include "command_line_parser.hpp"

#include "program_options.hpp"
#include "utils.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdio.h>
#include <string>
#include <sstream>
#include <sys/ioctl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <exception>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/any.hpp>
namespace po = boost::program_options;

Parser::Parser(int argc, char** argv)
{
    query_terminal_width();
    process_command_line(argc, argv);
}


void Parser::query_terminal_width()
{
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    terminal_width = 120 < w.ws_col ? 120 : w.ws_col;
}


Parser::~Parser()
{
}


void Parser::process_command_line(int argc, char** argv)
{
    // po::options_description desc("Allowed options. Please, read carefully for proper usage!!",
    //                              terminal_width, terminal_width/2);
    po::options_description desc("Allowed options. Please, read carefully for proper usage!!");
    try {
        desc.add_options()
        ("help", "Produce this help message.")
        ("skip-cbct",   po::bool_switch(&skip_cbct)->
                        default_value(skip_cbct),
                        "If the calculation on the CBCT should be skipped for debugging purposes.")
        // Common parameters
        ("plan_dir",    po::value<std::string>(&patient)->
                            required(),
                            "Base patient directory as produced by MCAuto-Astroid.")
        ("adapt_dir",   po::value<std::string>(&new_patient)->
                            required(),
                            "New plan directory.")
        ("cbct",        po::value<std::string>(&cbct_file)->
                            required(),
                            "CBCT in MHA format to adapt the plan to.")
        ("vf",          po::value<std::string>(&vf_file)->
                            required(),
                            "Vector field file from CT to CBCT. The file will be queried by "
                            "Plastimatch.")
        ("machine",     po::value<std::string>(&machine)->
                            default_value(machine),
                            "Machine to use during the gPMC simulation, beam model calculation or "
                            "adaptation energy limits. Case insensitive.")
        // Output files
        ("outplan",     po::value<std::string>(&out_plan)->
                            required(),
                            "Directory to write adapted gPMC inputs and final results.")
        ("optdir",      po::value<std::string>(&work_dir)->
                            required(),
                            "Directory to store optimization and plan evaluation files. "
                            "Will be prepended to any output if they don't contain \'/\'. This is "
                            "the working directory where subproducts will be stored.")
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
        ("vf_report",   po::value<std::string>(&vf_report_file)->
                            default_value(vf_report_file)->
                            implicit_value(implicit_vf_report_file),
                            "If a report should be generated with information about the vector "
                            "field and the subsequent spots movement. Requires out_vf and "
                            "out_shifts options.")
        // Weight adjustment modes
        ("method",      po::value<std::string>(&adapt_method_str)->
                            default_value(adapt_method_geometric_str),
                            "Adaptation method to use. The code is able to adapt a plan to a given "
                            "fraction using several methods. The basic approach consists of a the "
                            "geometrical adaptation focusing only in position and energy shifts, "
                            " without adapting spot weights. Additionally, it may calculate Dij "
                            "matrices with a native (naive) beam model or gPMC. Lastly, it can "
                            "also try to correct cold/hot spots after geometrical adaptation."
                            "Choices:\n"
                            "  -\"geometrical\": The spots are moved according to the passed "
                            "vector field. The energy is corrected according to is as well after "
                            "doing WEPL calculations. The vector field is probed at the end of "
                            "range position along the central axis of each spot. Although this may "
                            "not be exact in a patient geometry, the vector field is has been "
                            "observed to be smooth enough to allow for these uncertainty to be "
                            "negligible. All other methods are based on this as first step!!\n"
                            "  -\"beam_model\": If the Dij should be calculated with the native "
                            "(naive) beam model. The plan dose distribution is calculated at set "
                            "of positions within a mask at the plan and fraction geometry. It is "
                            "done in both for consistency purposes. A Dij is constructed at "
                            "fraction geometry and the weights are optimized to match the plan "
                            "dose distribution.\n"
                            "  -\"gpmc_dij\": A Dij matrix is calculated in a masked region with "
                            "gPMC. The plan dose distribution is read from the passed file and the "
                            "weights are optimized to match it\n"
                            "  -\"cold_spots\": gPMC is called on the geometrically adaptated "
                            "plan and the cold/hot areas are iteratively corrected.")
        // Weight adjustment files
        ("sf_dij",      po::value<float>(&spot_factor_dij),
                            "Spot factor for Dij calculation wiht gPMC.")
        ("sf_dose",     po::value<float>(&spot_factor_dose),
                            "Spot factor for dose calculation wiht gPMC.")
        ("plan_dose",   po::value<std::string>(&dose_plan_file),
                            "Plan dose distribution, usually from gPMC.")
        ("field_plan_doses", po::value<std::vector<std::string>>()->
                            multitoken(),
                            "Plan dose distributions per field warped to be used as ground truth, "
                            "Required if using cold spots adaptation. As output by gPMC or in mha.")
        ("frac_dose",   po::value<std::string>(&dose_frac_file)->
                            default_value(default_dose_frac_file)->
                            implicit_value(default_dose_frac_file),
                            "File root of fraction dose distribution after geometrical adaptation. "
                            "\".beam_#.dose\" will be appended to the string here passed. This is "
                            "an output file calculated with gPMC and then read for cold/hot spot "
                            "corrections.")
        ("plan_dij",    po::value<std::string>(&dij_plan_file),
                            "Dij file corresponding to the initial plan. Can be interpreted as "
                            "an output file in beam_model mode or as an input file. "
                            "Required when using beam model.")
        ("frac_dij",    po::value<std::string>(&dij_frac_file)->
                            default_value(default_dij_frac_file)->
                            implicit_value(default_dij_frac_file),
                            "Dij file name root corresponding to the daily geometry plan. In the "
                            "\"gpmc_dij\" it will be modified by appending \"beam_#.dij\". It is "
                            "calculated in the daily geometry and after applying the vector field "
                            "to the initial plan. Can be an output from the \"beam_model\" or the "
                            "\"gpmc_dij\" mode. Not necessary for geometric mode.")
        ("vf_masks",     po::value<std::vector<std::string>>()->
                            multitoken()->
                            implicit_value(std::vector<std::string>{""}, ""),
                            "List of binary images of a masks in CT space. They represent the "
                            "volume in which the vector field is probed. The order in which the "
                            "masks are passed represents their importance. Spots are raytraced in "
                            "CT space and the vector field is probed at the distal position along "
                            "the beamlet central axis (endpoint). If the beamlet enters any of the "
                            "masked regions, the enpoint will be defined as the last position "
                            "along the ray in the mask. If it doesn't, the endpoint will simply be "
                            "the last position along the central axis. If the ray enters a mask of "
                            "higher importance, then the endpoint will be defined w.r.t. that "
                            "region and other masks will be ignored. A usual example would be to "
                            "pass the target volume and 0.5, 1.0, 1.5 and 2.0 cm rings around it. "
                            "If they ray enters the target volume then the endpoint will be "
                            "w.r.t. to it. If it doesn't, but it enters the 0.5 cm, then that will "
                            "be the reference region. This all is important because the spots are "
                            "moved according to the vector field value at the endpoints. If the "
                            "command is not provided, the vector field will not be masked. An "
                            "additional explicit syntaxes to set the importance levels are given "
                            "at examples 2, 3, 4 and 5. 1 and 2 and 4 and 5 are equivalent.\n"
                            "Example:\n"
                            "   1) target.mha target+1cm.mha target+2cm.mha\n"
                            "   2) 1 target.mha 2 target+1cm.mha 3 target+2cm.mha\n"
                            "   3) 1 target.mha 1 target+1cm.mha 2 target+2cm.mha\n"
                            "   4) target.mha target+1cm.mha target+2cm.mha 3 target+3cm.mha\n"
                            "   5) 1 target.mha 2 target+1cm.mha 3 target+2cm.mha 3 target+3cm.mha")
        ("skin_mask",   po::value<std::vector<std::string>>()->
                            multitoken()->
                            implicit_value(std::vector<std::string>{""}, ""),
                            "Binary image of one or more masks in CT to be used during VF probing. "
                            "It is used as secondary criteria to limit the VF probing volume. If a "
                            "ray has been trace through the target mask passed to \"vf_mask\", "
                            "this is mask is ignored.")
        ("target_mask", po::value<std::vector<std::string>>()->
                            multitoken(),
                            "Binary image of one or more masks in CBCT to be used during Dij "
                            "calculation with the beam model, to feed to gPMC for Dij calculation, "
                            "to evaluate the dose on or to feed to Opt4D as target for "
                            "reoptimization.")
        ("target_rim",  po::value<std::vector<std::string>>()->
                            multitoken(),
                            "Binary image of one or more mask in CBCT. It is used by Opt4D to try to "
                            "keep the dose low outside the target.")
        ("oars",        po::value<std::vector<std::string>>()->
                            multitoken(),
                            "Binary image of one or more mask of OARs in CBCT. Used by Opt4D to keep "
                            "the dose low in these regions.")
        ("presciption", po::value<float>(&dose_prescription)->
                            default_value(dose_prescription),
                            "Dose prescribed to the target. The current implementation only allows "
                            "a single dose level and will be applied at the composite area covered "
                            "by to the \"mask\".")
        // Launchers
        ("opt4D",       po::bool_switch(&launch_opt4D)->
                            default_value(launch_opt4D),
                            "If Opt4D should be launched. If a file destination is not provided as "
                            "well, the default will be set.")
        ("sim_adapt",   po::bool_switch(&launch_adapt_simulation)->
                            default_value(launch_adapt_simulation),
                            "If gPMC should be launched to simulate the adaptated plan.")
        // Adaptation methods
        ("constraint",  po::value<std::vector<std::string>>()->
                            multitoken()->
                            implicit_value(std::vector<std::string>{adapt_constraint_free_str},
                                           adapt_constraint_free_str),
                            "Constrains set on the adaptation process. One or more can be "
                            "activated. In case conflicting options are found, the most "
                            "restrictive ones will be selected. Multitoken option."
                            "Choices:\n"
                            "  -\"free\": Default choice. Spots are moved freely in position and "
                            "energy.\n"
                            "  -\"fixed\": Spot positions and energies are not adjusted to the new "
                            "geometry.\n"
                            "  -\"fixed_pos\": Spots are not moved but their energy is adjusted to "
                            "the new geometry.\n"
                            "  -\"range_shifter\": A field-specific range shifter is calculated as "
                            "the average energy shift of the deepest layers in the field. The "
                            "deepest layers are defined as those that have at least the average "
                            "spot energy of the field. If the calculated range shifter has "
                            "positive thickness, a physical range shifter is imposed in the gPMC "
                            "simulation. In case the plan has a range shifter, it will be adapted. "
                            "The range shifter adaptation is done by changing it's density in the "
                            "simulation. This is an approximation done to keep the position of the "
                            "devices constant and avoid overlaps. Any range shift that can not be "
                            "corrected by modifying the density of the physical device will fall "
                            "back to the modification of the spot energies to reflecting any "
                            "remaining range shift, that is, a virtual range shifter is used.\n"
                            "  -\"virt_range_shifter\": A virtual range shifter is added following "
                            "the same procedure explained for the range shifter. In this case, the "
                            "corresponding energy is added or substracted from the spots "
                            "themselves, instead of modifying a physical device.\n"
                            "  -\"iso_shift\": A global isocenter shift is calculated as the "
                            "average displacement.\n"
                            "  -\"iso_shift_field\": A field-specific isocenter shift.")
        ("rshifter_steps", po::value<std::vector<std::string>>()->
                            multitoken()->
                            default_value(std::vector<std::string>{rshifter_steps_half_cm},
                                          rshifter_steps_half_cm)->
                            implicit_value(std::vector<std::string>{rshifter_steps_half_cm},
                                           rshifter_steps_half_cm),
                            "Possible WEPL values taken by modified or new range shifters. "
                            "Thicknesses are discretized by default. This is common practice in "
                            "the clinic because range shifters are commissioned and there are "
                            "limited devices. On the other hand, the range to energy conversion in "
                            "the range shifter calculator algorithm is a non-invertible fit to the "
                            "range of topasmediumspots beamlets in water. This means that there is "
                            "another fit to convert from energy to range that will yield slightly "
                            "different values and may cause the range shifter to calculate its "
                            "thickness eternally when the energy of a spot is outside the gantry "
                            "energy range. This is solved by the discretization. The possible "
                            "options are a set of predefined thicknesses, a list of WEPLs given by "
                            "the user or an interval as explained below."
                            "Choices:\n"
                            "  -\"mgh\": 3.45, 5.405 and 9.2 cm of water, corresponding to 0, 3, "
                            "4.7, 8 cm of PMMA 1.15 g/cm3.\n"
                            "  -\"cm\": 1, 2, 3 ... cm of water.\n"
                            "  -\"half_cm\": 0.5, 1.0, 1.5, 2.0 ... cm of water.\n"
                            "  -\"free\": No discretization.\n"
                            "  -List of values : \"1.234 2.345 3.456\" ... (space-separated).\n"
                            "  -Interval: \"x0.75\" will produce \"0.75 1.50 2.25\" ... .\n"
                            );

        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (argc == 1)
            throw std::invalid_argument("No input arguments.");
        if (vm.count("help") || argc == 1) {
            std::cerr << desc << std::endl;
            exit(EXIT_SUCCESS);
        }
        po::notify(vm);


        // MULTITOKEN OPTIONS
        if (!vm["target_mask"].empty())
            target_mask_files = vm["target_mask"].as<std::vector<std::string>>();
        if (!vm["target_rim"].empty())
            target_rim_files = vm["target_rim"].as<std::vector<std::string>>();
        if (!vm["oars"].empty())
            oars_files = vm["oars"].as<std::vector<std::string>>();
        std::vector<std::string> masks_command;
        if (!vm["vf_masks"].empty())
            masks_command = vm["vf_masks"].as<std::vector<std::string>>();
        if (!vm["skin_mask"].empty())
            skin_mask_files = vm["skin_mask"].as<std::vector<std::string>>();
        if (!vm["constraint"].empty()) {
            constraint_vec = vm["constraint"].as<std::vector<std::string>>();
            for (size_t i=0; i<constraint_vec.size(); ++i) {
                std::string str = constraint_vec.at(i);
                auto pos = str.find(' ');
                if (pos != std::string::npos) {
                    // replace multiple spaces with single
                    std::regex re("  +");
                    str = std::regex_replace(str, re, " ");
                    // remove leading and trailing spaces
                    if (pos == 0)
                        str = str.substr(1, str.size()-1);
                    if (str.find_last_of(' ') == str.size()-1)
                        str = str.substr(0, str.size()-2);
                    // append new values to vector
                    std::istringstream iss(str);
                    for(std::string s; iss >> s;)
                        constraint_vec.push_back(s);
                    // remove old element with spaces
                    constraint_vec.erase(constraint_vec.begin()+i);
                }
            }
        }
        if (!vm["rshifter_steps"].empty())
            rshifter_steps_vec = vm["rshifter_steps"].as<std::vector<std::string>>();
        if (!vm["field_plan_doses"].empty())
            field_dose_plan_files = vm["field_plan_doses"].as<std::vector<std::string>>();

        // PARSE INPUT MASKS FOR VF
        if (masks_command.size()) {
            std::regex rgx("[0-9]+");
            for (size_t i = 0; i < masks_command.size(); ++i) {
                if (std::regex_match(masks_command.at(i), rgx)) {
                    v_field_mask_importances.push_back(std::stoi(masks_command.at(i)));
                } else {
                    v_field_mask_files.push_back(masks_command.at(i));
                    if (v_field_mask_importances.size() < v_field_mask_files.size()) {
                        if (i == 0) {
                            v_field_mask_importances.push_back(1);
                        } else {
                            v_field_mask_importances.push_back(v_field_mask_importances.back()+1);
                        }
                    }
                }
            }
        }

        // ADAPT METHOD OPTIONS
        adapt_method_str = utils::toLower(adapt_method_str);
        if (adapt_method_str == adapt_method_geometric_str) {
            adapt_method = Adapt_methods_t::GEOMETRIC;
        } else if (adapt_method_str == adapt_method_beam_model_str) {
            adapt_method = Adapt_methods_t::BEAM_MODEL;
            if (dij_frac_file.empty()) {
                throw std::invalid_argument("Name of output fraction Dij file is necessary in "
                                            "\"beam model\" mode.");
            } else if (dij_plan_file.empty()) {
                std::cerr << "WARNING! No output given as Dij at plan stage, "
                             "required in beam model mode. " << std::endl;
                std::cerr << "Deduced from Dij at fraction as: ";
                std::string ext = "." + utils::get_file_extension(dij_frac_file);
                dij_plan_file = utils::replace_string(dij_frac_file, ext, "_at_plan"+ext);
                std::cerr << dij_plan_file << std::endl;
            } else if (target_mask_files.size() == 0) {
                throw std::invalid_argument("No mask given to reduce estimation space of the Dij.");
            }
        } else if (adapt_method_str == adapt_method_gpmc_dij_str) {
            adapt_method = Adapt_methods_t::GPMC_DIJ;
            if (dose_plan_file.empty()) {
                throw std::invalid_argument("the plan dose is required when using gPMC for "
                                            "Dij reoptimization");
            } else if (dij_frac_file.empty()) {
                throw std::invalid_argument("Name of output fraction Dij file is required when "
                                            "using gPMC for Dij reoptimization");
            }  else if (target_mask_files.size() == 0) {
                throw std::invalid_argument("No mask given to reduce calculation space of the Dij.");
            }
        } else if(adapt_method_str == adapt_method_cold_spots_str) {
            adapt_method = Adapt_methods_t::GPMC_DOSE;
            if (target_mask_files.size() == 0) {
                throw std::invalid_argument("No target mask to perform cold/hot spots correction.");
            } else if (target_rim_files.size() == 0) {
                throw std::invalid_argument("No target rim to reduce dose outside the target.");
            } else if (field_dose_plan_files.size() == 0) {
                throw std::invalid_argument("No plan dose per field were provided!");
            }
        } else {
            throw std::invalid_argument("Unrecognized adaptation method \""+adapt_method_str+"\"!");
        }
       
        // ADAPT METHOD CONSTRAINTS
        bool free = true;
        bool fixed = false;
        bool fixed_pos = false;
        bool v_range_shifter = false;
        bool range_shifter = false;
        bool iso_shift = false;
        bool iso_shift_field = false;
        for (auto str: constraint_vec) {
            str = utils::toLower(str);
            if (str == adapt_constraint_free_str) {
                free = true;
            } else if (str == adapt_constraint_fixed_str) {
                fixed = true;
            } else if (str == adapt_constraint_fixed_pos_str) {
                fixed_pos = true;
            } else if (str == adapt_constraint_rs_str) {
                range_shifter = true;
            } else if (str == adapt_constraint_v_rs_str) {
                v_range_shifter = true;
            } else if (str == adapt_constraint_iso_str) {
                iso_shift = true;
            } else if (str == adapt_constraint_iso_field_str) {
                iso_shift_field = true;
            } else {
                throw std::invalid_argument("Unrecognized adaptation constraint \""+str+"\"!");
            }
        }
        if (iso_shift && range_shifter)
            adapt_constraints = Adapt_constraints_t::ISOCENTER_SHIFT_RANGE_SHIFTER;
        else if (iso_shift_field && range_shifter)
            adapt_constraints = Adapt_constraints_t::FIELD_ISOCENTER_SHIFT_RANGE_SHIFTER;
        else if (iso_shift)
            adapt_constraints = Adapt_constraints_t::ISOCENTER_SHIFT;
        else if (iso_shift_field)
            adapt_constraints = Adapt_constraints_t::FIELD_ISOCENTER_SHIFT;
        else if (range_shifter)
            adapt_constraints = Adapt_constraints_t::RANGE_SHIFTER;
        else if (iso_shift && v_range_shifter)
            adapt_constraints = Adapt_constraints_t::ISOCENTER_SHIFT_V_RANGE_SHIFTER;
        else if (iso_shift_field && v_range_shifter)
            adapt_constraints = Adapt_constraints_t::FIELD_ISOCENTER_SHIFT_V_RANGE_SHIFTER;
        else if (v_range_shifter)
            adapt_constraints = Adapt_constraints_t::V_RANGE_SHIFTER;
        else if (fixed_pos)
            adapt_constraints = Adapt_constraints_t::FIXED_POS;
        else if (fixed)
            adapt_constraints = Adapt_constraints_t::FIXED;
        else if (free)
            adapt_constraints = Adapt_constraints_t::FREE;

        // RANGE SHIFTER STEPS
        if (rshifter_steps_vec.size() == 1) {
            std::string temp = utils::toLower(rshifter_steps_vec.front());
            if (temp == rshifter_steps_free) {
                rshifter_steps.set_mode(RShifter_steps_t::FREE);
            } else if (temp == rshifter_steps_half_cm) {
                rshifter_steps.set_mode(RShifter_steps_t::HALF_CM);
            } else if (temp == rshifter_steps_cm) {
                rshifter_steps.set_mode(RShifter_steps_t::CM);
            } else if (temp == rshifter_steps_mgh) {
                rshifter_steps.set_mode(RShifter_steps_t::MGH);
            } else if (temp.at(0) == 'x') {
                rshifter_steps.set_mode(RShifter_steps_t::INTERVAL, rshifter_steps_vec);
            } else {
                throw std::invalid_argument("Unrecognized range shifter discretization mode. It "
                                            "could not be deduced from \""+temp+"\"");
            }
        } else if (rshifter_steps_vec.size() > 1) {
            rshifter_steps.set_mode(RShifter_steps_t::LIST, rshifter_steps_vec);
        } else {
            std::string temp;
            for (auto i: rshifter_steps_vec)
                temp += i+" ";
            throw std::invalid_argument("Unrecognized range shifter discretization mode. It "
                                        "could not be deduced from \""+temp+"\"");
        }

        // PARSE PASSED TARGET MASKS AND VF AND VERIFY THEY EXIST
        target_expanded_files.reserve(target_mask_files.size() + target_rim_files.size());
        for (auto& str: target_mask_files) {
            if (std::ifstream(str)) {
                str = utils::get_full_path(str);
                target_expanded_files.push_back(str);
            } else {
                std::cerr << "Mask file not accessible: \""+str+"\"!" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        for (auto& str: target_rim_files) {
            if (std::ifstream(str)) {
                str = utils::get_full_path(str);
                target_expanded_files.push_back(str);
            } else {
                std::cerr << "Mask file not accessible: \""+str+"\"!" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        for (auto& str: oars_files) {
            if (std::ifstream(str)) {
                str = utils::get_full_path(str);
            } else {
                std::cerr << "Mask file not accessible: \""+str+"\"!" << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        // PARSE VF MASKS
        if (std::ifstream(vf_file)) {
            vf_file = utils::get_full_path(vf_file);
        } else {
            std::cerr << "VF file not accessible: \""+vf_file+"\"!" << std::endl;
        }
        if (v_field_mask_files.size() != 0) {
            for (auto& str: v_field_mask_files) {
                if (std::ifstream(str)) {
                    str = utils::get_full_path(str);
                } else {
                    std::cerr << "VF mask file not accessible: \""+str+"\"!" << std::endl;
                }
            }
        }
        if (skin_mask_files.size() != 0) {
            for (auto& str: skin_mask_files) {
                if (std::ifstream(str)) {
                    str = utils::get_full_path(str);
                } else {
                    std::cerr << "Skin mask file not accessible: \""+str+"\"!" << std::endl;
                }
            }
        }

        // REPORTING
        if ( !vf_report_file.empty() && ( data_shifts_file.empty() || data_vf_file.empty()))
            throw std::invalid_argument("report option needs shifts and vf output");

        // CORRECT PATHS
        new_patient = utils::get_full_path(new_patient);
        if (!work_dir.empty() &&
            !data_shifts_file.empty() && data_shifts_file.find('/') == std::string::npos)
            data_shifts_file = work_dir + '/' + data_shifts_file;
        if (!work_dir.empty() &&
            !data_vf_file.empty() && data_vf_file.find('/') == std::string::npos)
            data_vf_file = work_dir + '/' + data_vf_file;
        if (!work_dir.empty() &&
            !ct_traces_file.empty() && ct_traces_file.find('/') == std::string::npos)
            ct_traces_file = work_dir + '/' + ct_traces_file;
        if (!work_dir.empty() &&
            !cbct_traces_file.empty() && cbct_traces_file.find('/') == std::string::npos)
            cbct_traces_file = work_dir + '/' + cbct_traces_file;

        // CREATE OUTPUT DIRECTORIES
        mkdir(out_plan.c_str(), 0774);
        mkdir(work_dir.c_str(), 0774);
        
        // NORMALIZE INPUTS
        machine = utils::toLower(machine);
    } catch(std::exception& e) {
        std::cerr << "ERROR! " << e.what() << std::endl;
        print_inputs();
        std::cerr << desc << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Parser::print_inputs () {
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
        } else { // Assumes that the only remainder type is vector<string>
            try {
                std::vector<std::string> vect = vm[it->first].as<std::vector<std::string>>();
                uint i = 0;
                for (std::vector<std::string>::iterator oit=vect.begin();
                     oit != vect.end(); oit++, ++i) {
                    std::cout << "\r> " << it->first << "[" << i << "] = " << (*oit) << std::endl;
                }
            } catch (const boost::bad_any_cast &) {
                std::cout << "UnknownType(" << ((boost::any)it->second.value()).type().name();
                std::cout << ")" << std::endl;
            }
        }
    }
    std::cout << std::endl;
}
