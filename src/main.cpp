#include <iostream>
#include <iomanip>
#include <string>
#include <limits>
#include <boost/timer/timer.hpp>

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
    // Timer -----------------------------------------------------------------
    boost::timer::cpu_times time_ct_tracing, time_cbct_tracing;
    boost::timer::cpu_times time_geometric_sim, time_opt4d, time_opt4d_validation;
    boost::timer::cpu_timer timer;

    // Read input parameters -------------------------------------------------
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

    // 1: Load CT
    Volume_t ct(pat.planning_ct_file, Volume_t::Source_type::CTVOLUME);
    ct.setDims(pat.ct); // The CT volume lacks dimensions information
    // 2: Get endpoints and warped endpoints in CT
    compute_in_ct (pat, parser, ct,                   // inputs
                   warper,                            // outputs
                   endpoints, initpos,                // outputs
                   warped_endpoints, warped_initpos); // outputs
    // 4: Unload CT (I could destroy the object, but that could raise problems)
    ct.freeMemory();

    time_ct_tracing = timer.elapsed();
    timer.start();

    // Stop process for debug purposes
    if (parser.skip_cbct)
        exit(EXIT_SUCCESS);

    // 5: Load CBCT
    if (parser.adapt_constraints != Adapt_constraints_t::FIXED) {
        Volume_t cbct(parser.cbct_file, Volume_t::Source_type::MHA);
        // pat.update_geometry(cbct);
        // 6: Get endpoints and corrections in CBCT
        Array3<float> traces_data;
        compute_in_cbct (pat, parser, cbct,                // inputs
                         warped_endpoints, warped_initpos, // inputs
                         new_energies, cbct_endpoints, traces_data); // outputs
        // 8: Unload CBCT (I could destroy the object, but that could raise problems)
        cbct.freeMemory();
        // 7: constrain energy shifts and export geometric adaptation (weight_scaling = [1 .. 1])
        calculate_range_shifter(parser.rshifter_steps, parser.adapt_constraints, new_energies,
                                pat.range_shifters, pat.machine, pat.isocenter_to_beam_distance,
                                pat.source_energies, pat.spots_per_field);
    } else {
        new_energies = pat.source_energies;
        std::transform(new_energies.begin(), new_energies.end(), new_energies.begin(),
            std::bind1st(std::multiplies<float>(), 1e6));
    }

    time_cbct_tracing = timer.elapsed();
    timer.start();

    export_adapted (pat, parser.work_dir, parser.data_shifts_file,
                    new_energies, weight_scaling,
                    warped_initpos, warped_endpoints, warper,
                    pat.geometric_tramp_names);

    // 9: Free memory in device
    freeCTMemory();

    if (!parser.vf_report_file.empty()) {
        generate_report(parser.vf_report_file, parser.data_vf_file, parser.data_shifts_file,
                        parser.out_plan, pat.tramp_files);
    }

    // 10: Output files for gPMC simulation
    if (parser.adapt_method == Adapt_methods_t::GEOMETRIC) {
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

        time_geometric_sim = timer.elapsed();
        timer.start();

    } else if (parser.adapt_method == Adapt_methods_t::GPMC_DOSE) {
        Gpmc_manager gpmc_dij(pat, parser.new_patient, parser.dij_frac_file, "dosedij",
                              parser.out_plan, parser.work_dir, pat.adapted_tramp_names, warper.vf_ave);
        gpmc_dij.write_dij_files(parser.spot_factor_dij, 0, true,
                                 utils::join_vectors(parser.target_mask_files, parser.oars_files,
                                                     parser.target_rim_files));
        // gpmc_dij.launch();

        time_geometric_sim = timer.elapsed();
        timer.start();

        adapt_weights(gpmc_dij.get_field_dij_files(),
                      parser.target_mask_files, parser.target_rim_files, parser.oars_files,
                      parser.dose_prescription, gpmc_dij.get_to_Gy_factor(),
                      parser.work_dir,
                      pat.spots_per_field,
                      pat.source_weights, weight_scaling);
        time_opt4d = timer.elapsed();
        timer.start();

        if (!parser.vf_report_file.empty()) {
            generate_report(parser.vf_report_file, parser.data_vf_file, parser.data_shifts_file,
                            parser.out_plan, pat.tramp_files);
        }

        // 11: Export results and report
        export_adapted (pat, parser.out_plan, parser.data_shifts_file,
                        new_energies, weight_scaling,
                        warped_initpos, warped_endpoints, warper,
                        pat.adapted_tramp_names);

        // Check adaptation!!
        Gpmc_manager gpmc_dose(pat, parser.new_patient, parser.dose_frac_file, "dose",
                               parser.out_plan, parser.out_plan, pat.adapted_tramp_names,
                               warper.vf_ave);
        gpmc_dose.write_dose_files(parser.spot_factor_dose);
        if (parser.launch_adapt_simulation) {
            gpmc_dose.launch();
            time_opt4d_validation = timer.elapsed();
            timer.start();
        }
    }

    // Report time ------------------------------------------------------------
    std::cout << std::endl;
    std::cout << "Total time:" << std::endl;
    std::cout << "    CT tracing:     " << time_ct_tracing.wall << std::endl;
    std::cout << "    CBCT tracing:   " << time_cbct_tracing.wall << std::endl;
    std::cout << "    gPMC geometric: " << time_geometric_sim.wall << std::endl;
    std::cout << "    Opt4D:          " << time_opt4d.wall << std::endl;
    std::cout << "    Opt validation: " << time_opt4d_validation.wall << std::endl;

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
    // I use std::stringstream to perform only one write (<<) to std::cout. It is guaranteed to 
    // be atomic in C++11, but that applies to each << call. By using stringstreams, only one << 
    // is required. I finally did not use OpenMP, but I'll still leave this implementation
    std::stringstream stream;
    stream << "Generating report ..." << std::endl;
    std::string interp = "python3 ";
    std::string code   = std::string(INSTALLATION_PATH) + "/src/extra/create_adapt_report.py ";
    std::string vf     = "--vf " + data_vf_file + " ";
    std::string shifts = "--shifts " + data_shifts_file + " ";
    std::string tramps = "--tramps ";
    for (size_t i = 0; i < tramp_files.size(); i++)
        tramps += tramp_files.at(i) + " ";
    std::string outdir = "--outdir " + out_dir + " ";
    std::string outfile = "--outfile " + vf_report_file;
    std::string command = interp + code + vf + shifts + tramps + outdir + outfile + " &";
    stream << "Running: " << command << std::endl;
    std::cout << stream.str();
    int res = system(command.c_str());

    if (res) {
        stream.str(std::string());
        stream.clear();
        stream << "ERROR! Reporting failed with error code: " << res << std::endl;
        std::cerr << stream.str();
        exit(EXIT_FAILURE);
    }
}
