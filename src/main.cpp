#include <iostream>
#include <iomanip>
#include <string>
#include <limits>

#include "command_line_parser.hpp"
#include "cold_spots_fixer.hpp"
#include "influence_manager.hpp"
#include "initialize_rays.cuh"
#include "enviroment.hpp"
#include "gpmc_manager.hpp"
#include "gpu_ct_to_device.cuh"
#include "gpu_main.cuh"
#include "gpu_source_positioning.cuh"
#include "opt4D_manager.hpp"
#include "patient_parameters.hpp"
#include "program_options.hpp"
#include "tramp.hpp"
#include "structure_sampler.hpp"
#include "utils.hpp"
#include "volume.hpp"
#include "warper.hpp"

// DONE Read CT
// DONE Raytrace endpoints on CT
// DONE Apply VF to endpoints -> updated endpoints
// DONE Apply VF to starting points -> updated starting points
// DONE Projection onto treatment plane
// DONE Read CBCT
// DONE Raytrace vf_start_points in CBCT
// DONE Get energy distance to updated endpoints
// DONE Convert data back to tramp file
// DONE Fix mha output's offset to visualize in Slicer
// DONE Verify that VF probing is correctly done. How?
// Done Verify CBCT offsets are correctly updated

void compute_in_ct(const Patient_Parameters_t& pat,
                   const Parser& parser,
                   const Volume_t& ct,
                   Warper_t& warper,
                   Array4<float>& endpoints,
                   Array4<float>& initpos,
                   Array4<float>& warped_endpoints,
                   Array4<float>& warped_initpos);

void compute_in_cbct(Patient_Parameters_t& pat,
                     const Parser& parser,
                     const Volume_t& cbct,
                     const Array4<float>& ct_warped_endpoints,
                     const Array4<float>& ct_warped_initpos,
                     std::vector<float>& new_energies,
                     Array4<float>& cbct_endpoints,
                     Array3<float>& traces_data);

// void process_influences (Array4<float> influence_ct,
//                          Array4<float> influence_cbct);

void export_adapted(Patient_Parameters_t& pat,
                    const std::string& out_dir,
                    const std::string& shifts_file,
                    const std::vector<float>& energy_shift,
                    const std::vector<float>& weight_scaling,
                    Array4<float> pat_pos,
                    Array4<float> pat_pos2,
                    const Warper_t& warper,
                    std::vector<std::string> tramp_names = std::vector<std::string>());

void export_shifts(const std::vector<float> e,
                   const std::vector<float>& w,
                   const Array4<float>& p,
                   const std::string& file,
                   const short& beamid,
                   const Vector3_t<float>& isocenter_shift);

void generate_report(const std::string& vf_report_file,
                     const std::string& data_vf_file,
                     const std::string& data_shifts_file,
                     const std::string& out_dir,
                     const std::vector<std::string>& tramp_files);

int main(int argc, char** argv)
{
    Parser parser(argc, argv);
    parser.print_inputs();

    // Read input patient ----------------------------------------------------
    Patient_Parameters_t pat(parser.patient);
    pat.print();
    pat.ct_to_int_coordinates();
    pat.set_treatment_planes();

    // Start device ----------------------------------------------------------
    cudaEvent_t start;
    initialize_device(start);

    // Data containers -------------------------------------------------------
    Array4<float> endpoints(pat.total_spots);
    Array4<float> initpos(pat.total_spots);
    Array4<float> warped_endpoints(pat.total_spots);
    Array4<float> warped_initpos(pat.total_spots);
    Array4<float> cbct_endpoints(pat.total_spots);
    std::vector<float> new_energies(pat.total_spots);
    std::vector<float> weight_scaling(pat.total_spots, 1.0);

    // Warper ----------------------------------------------------------------
    Warper_t warper(parser.vf_file, parser.data_vf_file);
    warper.print_vf();

    // Adaptive steps --------------------------------------------------------
    // Initialize influence manager


    // Array4<float> influence_ct(pat.total_spots*pat.total_spots);
    // Array4<float> influence_cbct(pat.total_spots*pat.total_spots);
    // uint nprobes = std::numeric_limits<uint>::max();
    // uint nprobes = pat.total_spots;
    // structure_sampler (parser.ct_mask_file, percentage, pat.total_spots, warper, pat.ct,
                       // influence_ct, influence_cbct);

    // 1: Load CT
    Volume_t ct(pat.planning_ct_file, Volume_t::Source_type::CTVOLUME);
    ct.setDims(pat.ct); // The CT volume lacks dimensions information
    // 2: Get endpoints and warped endpoints in CT
    compute_in_ct (pat, parser, ct,                   // inputs
                   warper,                            // outputs
                   endpoints, initpos,                // outputs
                   warped_endpoints, warped_initpos); // outputs
    // 3: Get influence in CT: I need the vf average and that is calculated in compute_in_ct
    // Influence_manager influences(parser, pat, warper, ct.getMetadata());
    // influences.get_dose_at_plan();
    // 4: Unload CT (I could destroy the object, but that could raise problems)
    ct.freeMemory();

    // Stop process for debug purposes
    if (parser.skip_cbct)
        exit(EXIT_SUCCESS);

    // 5: Load CBCT
    Volume_t cbct(parser.cbct_file, Volume_t::Source_type::MHA);
    // pat.update_geometry(cbct);
    // 6: Get endpoints and corrections in CBCT
    Array3<float> traces_data;
    compute_in_cbct (pat, parser, cbct,                // inputs
                     warped_endpoints, warped_initpos, // inputs
                     new_energies, cbct_endpoints, traces_data); // outputs
    // 7: constrain energy shifts and export geometric adaptation (weight_scaling = [1 .. 1])
    calculate_range_shifter(parser.rshifter_steps, parser.adapt_constraints, new_energies,
                            pat.range_shifters, pat.machine, pat.isocenter_to_beam_distance,
                            pat.source_energies, pat.spots_per_field);
    export_adapted (pat, parser.work_dir, parser.data_shifts_file,
                    new_energies, weight_scaling,
                    warped_initpos, warped_endpoints, warper,
                    pat.geometric_tramp_names);
    // 7: Get influence in CBCT
    // influences.get_dij_at_frac (new_energies);
    // 8: Unload CBCT (I could destroy the object, but that could raise problems)
    cbct.freeMemory();

    // 9: Free memory in device
    freeCTMemory();

    // 10: Correct weights
    // if (parser.launch_opt4D) {
    //     Opt4D_manager opt4d(parser.work_dir);
    //     opt4d.populate_directory(influences.n_spots, influences.n_voxels,
    //                              influences.matrix_at_plan, influences.matrix_at_frac);
    //     opt4d.launch_optimization();
    //     weight_scaling = opt4d.get_weight_scaling();
    // }

    // 10: Output files for gPMC simulation
    if (parser.adapt_method == Adapt_methods_t::GEOMETRIC){
        Gpmc_manager gpmc(pat, parser.new_patient, parser.dose_frac_file, "dose",
                          parser.out_plan, parser.work_dir, pat.adapted_tramp_names, warper.vf_ave);
        gpmc.write_dose_files(parser.spot_factor_dose);
        if (parser.launch_adapt_simulation) {
            gpmc.launch();
        }
        Volume_t target_mask = utils::read_masks (parser.target_mask_files);
        Volume_t oars_mask = utils::read_masks (parser.oars_files);
        check_adaptation_from_dose(gpmc.get_total_dose_file(), target_mask, oars_mask,
                                   parser.dose_prescription, gpmc.get_to_Gy_factor());
    } else if (parser.adapt_method == Adapt_methods_t::GPMC_DOSE) {
        Gpmc_manager gpmc_dij(pat, parser.new_patient, parser.dij_frac_file, "dosedij",
                              parser.out_plan, parser.work_dir, pat.adapted_tramp_names, warper.vf_ave);
        gpmc_dij.write_dij_files(parser.spot_factor_dij, 0, true,
                                 utils::join_vectors(parser.target_mask_files, parser.oars_files,
                                                     parser.target_rim_files));
        gpmc_dij.launch();

        std::valarray<float> adapt_dose;
        std::valarray<float> adapt_dose_in_mask;
        std::valarray<float> adapt_dose_in_target;
        std::vector<std::valarray<float>> adapt_field_dose;
        std::vector<std::valarray<float>> adapt_field_dij;
        std::valarray<float> underdose_mask;
        Volume_t target_mask = utils::read_masks (parser.target_mask_files);
        Volume_t target_rim_mask = utils::read_masks (parser.target_rim_files);
        Volume_t oars_mask = utils::read_masks (parser.oars_files);
        check_adaptation(gpmc_dij.get_field_dij_files(),
                         parser.dose_prescription, gpmc_dij.get_to_Gy_factor(),
                         pat.spots_per_field, adapt_field_dij,
                         adapt_dose, adapt_field_dose, underdose_mask,
                         target_mask, target_rim_mask, oars_mask,
                         adapt_dose_in_mask, adapt_dose_in_target);
        // output_debug_doses (parser.work_dir, underdose_mask, target_mask,
        //                     parser.field_dose_plan_files,
        //                     adapt_dose, adapt_field_dose,
        //                     adapt_dose_in_target,
        //                     parser.dose_prescription);
        cold_spots_fixer (adapt_dose_in_mask,
                          adapt_field_dij, target_mask,
                          target_rim_mask, oars_mask,
                          parser.dose_prescription,
                          parser.work_dir, pat.spots_per_field,
                          pat.source_weights,
                          weight_scaling);

        // Check adaptation!!
        Gpmc_manager gpmc_dose(pat, parser.new_patient, parser.dose_frac_file, "dose",
                               parser.out_plan, parser.out_plan, pat.adapted_tramp_names,
                               warper.vf_ave);
        gpmc_dose.write_dose_files(parser.spot_factor_dose);
        if (parser.launch_adapt_simulation) {
            gpmc_dose.launch();
        }
    }

    // 11: Export results and report
    export_adapted (pat, parser.out_plan, parser.data_shifts_file,
                    new_energies, weight_scaling,
                    warped_initpos, warped_endpoints, warper,
                    pat.adapted_tramp_names);

    if (!parser.vf_report_file.empty()) {
        generate_report(parser.vf_report_file, parser.data_vf_file, parser.data_shifts_file,
                        parser.out_plan, pat.tramp_files);
    }

    // Stop device
    stop_device(start);
    exit(EXIT_SUCCESS);
}


void compute_in_ct(const Patient_Parameters_t& pat,
                   const Parser& parser,
                   const Volume_t& ct,
                   Warper_t& warper,
                   Array4<float>& endpoints,
                   Array4<float>& initpos,
                   Array4<float>& warped_endpoints,
                   Array4<float>& warped_initpos)
{
    // Get endpoints in CT ----------------------------
    Array4<float> initpos_xbuffer_dbg(pat.total_spots);
    gpu_raytrace_original (pat, parser, ct, endpoints, initpos_xbuffer_dbg, initpos);

    // Warp endpoints in CT ---------------------------
    warped_endpoints.resize(pat.total_spots);
    warped_initpos.resize(pat.total_spots);
    std::copy (endpoints.begin(), endpoints.end(), warped_endpoints.begin());
    std::copy (initpos.begin(), initpos.end(), warped_initpos.begin());
    warper.apply_to_plan(warped_endpoints, warped_initpos, pat.ct, pat.treatment_planes,
                         pat.angles, pat.spots_per_field, parser.adapt_constraints);

    // Print results
    std::cout << "Warped patient positions and wepl (example):" << std::endl;
    size_t iters = pat.total_spots < 5 ? pat.total_spots : 5;
    for (size_t i = 0; i < iters; i++) {
        initpos.at(i).print_as_3D(" -> ");
        warped_initpos.at(i).print_as_3D();
    }
}


void compute_in_cbct(Patient_Parameters_t& pat,
                     const Parser& parser,
                     const Volume_t& cbct,
                     const Array4<float>& warped_endpoints,
                     const Array4<float>& warped_initpos,
                     std::vector<float>& new_energy,
                     Array4<float>& cbct_endpoints,
                     Array3<float>& traces_info)
{
    // Get endpoints in CBCT --------------------
    cbct_endpoints.resize(pat.total_spots);
    gpu_raytrace_warped (pat, parser, cbct, warped_endpoints,
                         warped_initpos, cbct_endpoints, new_energy,
                         traces_info);
}

void export_adapted(Patient_Parameters_t& pat,
                    const std::string& out_dir,
                    const std::string& shifts_file,
                    const std::vector<float>& energy_shift,
                    const std::vector<float>& weight_scaling,
                    Array4<float> warped_initpos,
                    Array4<float> warped_endpoints,
                    const Warper_t& warper,
                    std::vector<std::string> tramp_names)
{
    // Go to virtual source pos
    treatment_plane_to_virtual_src (warped_initpos, warped_endpoints, pat, warper.vf_ave);

    // From virtual source to iso
    virtual_src_to_iso_pos(warped_initpos, pat.virtualSAD);
    
    // Assign to new spotmap
    for (size_t i = 0; i < pat.nbeams; i++) {
        size_t offset_a = (i != 0) ? pat.accu_spots_per_field.at(i-1) : 0;
        size_t offset_b = pat.accu_spots_per_field.at(i);
        std::vector<float> subset_energies;
        std::vector<float> subset_weights_scaling;
        Array4<float> subset_pat_pos;

        utils::subset_vector<float>(subset_energies, energy_shift, offset_a, offset_b);
        utils::subset_vector<float>(subset_weights_scaling, weight_scaling, offset_a, offset_b);
        utils::subset_vector< Vector4_t<float>>(subset_pat_pos, warped_initpos, offset_a, offset_b);

        Tramp_t tramp(pat.tramp_files.at(i));
        tramp.set_new_energies(subset_energies);
        tramp.scale_weights(subset_weights_scaling);
        tramp.set_pos(subset_pat_pos);
        std::string tramp_name = pat.tramp_files.at(i);
        if (tramp_names.size() != 0)
            tramp_name = tramp_names.at(i);
        tramp.to_file(tramp_name, out_dir);

        if (!shifts_file.empty())
            export_shifts(tramp.get_last_energy_shift_eV(), subset_weights_scaling,
                          subset_pat_pos, shifts_file, i, warper.vf_ave);
    }
}

void export_shifts(const std::vector<float> e,
                   const std::vector<float>& w,
                   const Array4<float>& p,
                   const std::string& file,
                   const short& beamid,
                   const Vector3_t<float>& isocenter_shift)
{
    std::ofstream ofs;
    if (beamid != 0) {
        ofs.open (file, std::ios::out | std::ios::app);
    } else {
        ofs.open (file, std::ios::out);
    }

    utils::check_fs(ofs, file, "to write adaptation shifts.");

    if (beamid == 0) {
        std::cout << "Writting adaptation shifts to " << file << std::endl;
        std::string txt = "Isocenter shift: " +
                          std::to_string(isocenter_shift.x) + ", " + 
                          std::to_string(isocenter_shift.y) + ", " + 
                          std::to_string(isocenter_shift.z) + " cm";
        std::cout << txt << std::endl;
        ofs << "# " << txt << '\n';
        ofs << "# w e x y z d beamid spotid\n";
    }

    for (size_t spotid = 0; spotid < e.size(); spotid++) {
        float vx = p.at(spotid).x;
        float vy = p.at(spotid).y;
        ofs << w.at(spotid) << " " << e.at(spotid) << " ";
        ofs << p.at(spotid).x << " " << p.at(spotid).y << " ";
        ofs << p.at(spotid).z << " " << std::sqrt(vx*vx + vy*vy) << " ";
        ofs << beamid << " " << spotid << "\n";
    }
}

void generate_report(const std::string& vf_report_file,
                     const std::string& data_vf_file,
                     const std::string& data_shifts_file,
                     const std::string& out_dir,
                     const std::vector<std::string>& tramp_files)
{
    std::cout << "Generating report ..." << std::endl;
    std::string interp = "python3 ";
    std::string code   = std::string(INSTALLATION_PATH) + "/src/extra/create_adapt_report.py ";
    std::string vf     = "--vf " + data_vf_file + " ";
    std::string shifts = "--shifts " + data_shifts_file + " ";
    std::string tramps = "--tramps ";
    for (size_t i = 0; i < tramp_files.size(); i++)
        tramps += tramp_files.at(i) + " ";
    std::string outdir = "--outdir " + out_dir + " ";
    std::string outfile = "--outfile " + vf_report_file;
    std::string command = interp + code + vf + shifts + tramps + outdir + outfile;
    std::cout << "Running: " << command << std::endl;
    int res = system(command.c_str());

    if (res) {
        std::cerr << "ERROR! Reporting failed with error code: " << res << std::endl;
        exit(EXIT_FAILURE);
    }
}
