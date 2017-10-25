#include <iostream>
#include <iomanip>
#include <string>

#include "command_line_parser.hpp"
#include "initialize_rays.cuh"
#include "enviroment.hpp"
#include "gpu_main.cuh"
#include "gpu_source_positioning.cuh"
#include "opt4D_manager.hpp"
#include "patient_parameters.hpp"
#include "program_options.hpp"
#include "tramp.hpp"
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
// TODO Fix mha output's offset to visualize in Slicer
// TODO Verify that VF probing is correctly done. How?
// Done Verify CBCT offsets are correctly updated

void compute_in_ct(const Patient_Parameters_t& pat,
                   const Parser& parser,
                   Warper_t& warper,
                   Array4<double>& endpoints,
                   Array4<double>& initpos,
                   Array4<double>& warped_endpoints,
                   Array4<double>& warped_initpos,
                   Array4<double>& influence_ct);

void compute_in_cbct(Patient_Parameters_t& pat,
                     const Parser& parser,
                     const Array4<double>& ct_warped_endpoints,
                     const Array4<double>& ct_warped_initpos,
                     std::vector<double>& energy_shift,
                     Array4<double>& cbct_endpoints,
                     Array4<double>& influence_cbct);

// void process_influences (Array4<double> influence_ct,
//                          Array4<double> influence_cbct);

void export_adapted(Patient_Parameters_t& pat,
                    const Parser& pars,
                    const std::vector<double>& energy_shift,
                    Array4<double>& pat_pos,
                    Array4<double>& pat_pos2,
                    const Warper_t& warper);

void export_shifts(const std::vector<double>& e,
                   const Array4<double>& p,
                   const std::string& file,
                   const short& beamid,
                   const Vector3_t<double>& isocenter_shift);

void generate_report(const std::string& report,
                     const std::string& output_vf,
                     const std::string& output_shifts,
                     const std::vector<std::string>& tramp_files);

int main(int argc, char** argv)
{
    Parser parser(argc, argv);
    parser.print_parameters();

    // Read input patient ----------------------------------------------------
    Patient_Parameters_t pat(parser.patient);
    pat.print();
    pat.ct_to_int_coordinates();
    pat.set_treatment_planes();

    // Start device ----------------------------------------------------------
    cudaEvent_t start;
    initialize_device(start);

    // Data containers and warper --------------------------------------------
    Array4<double> endpoints(pat.total_spots);
    Array4<double> initpos(pat.total_spots);
    Array4<double> warped_endpoints(pat.total_spots);
    Array4<double> warped_initpos(pat.total_spots);
    Array4<double> cbct_endpoints(pat.total_spots);
    Array4<double> influence_ct(pat.total_spots*pat.total_spots);
    Array4<double> influence_cbct(pat.total_spots*pat.total_spots);
    std::vector<double> energy_shift(pat.total_spots);

    // Warper ----------------------------------------------------------------
    Warper_t warper(parser.vf_file, parser.output_vf);

    // Adaptive steps --------------------------------------------------------
    // 1: Get endpoints in CT
    compute_in_ct (pat, parser,                       // inputs
                   warper,                            // outputs
                   endpoints, initpos,                // outputs
                   warped_endpoints, warped_initpos,  // outputs
                   influence_ct);                     // outputs

    if (!parser.skip_cbct)
    {
        // 2: Get endpoints in CBCT and compare with CT to correct energy
        compute_in_cbct (pat, parser,                      // inputs
                         warped_endpoints, warped_initpos, // inputs
                         energy_shift, cbct_endpoints,     // outputs
                         influence_cbct);                  // outputs
    }
    // 3: Correct weights
    if (!parser.output_opt4D_files.empty())
    {
        Opt4D_manager opt4d(parser.output_opt4D_files);
        opt4d.populate_directory(influence_ct, influence_cbct);
        if (parser.launch_opt4D)
            opt4d.launch_optimization();
    }

    // Export results and report
    export_adapted (pat, parser, energy_shift, warped_initpos, warped_endpoints, warper);
    if (!parser.report.empty())
        generate_report(parser.report, parser.output_vf, parser.output_shifts, pat.tramp_files);

    // Stop device
    stop_device(start);
    exit(EXIT_SUCCESS);
}


// void process_influences (Array4<double> influence_ct,
//                          Array4<double> influence_cbct)
// {
// }

void compute_in_ct(const Patient_Parameters_t& pat,
                   const Parser& parser,
                   Warper_t& warper,
                   Array4<double>& endpoints,
                   Array4<double>& initpos,
                   Array4<double>& warped_endpoints,
                   Array4<double>& warped_initpos,
                   Array4<double>& influence)
{
    // Read CT and launch rays
    Volume_t ct(pat.planning_ct_file, Volume_t::Source_type::CTVOLUME);
    // The CT volume lacks dimensions information
    ct.setDims(pat.ct);

    // Get endpoints in CT ----------------------------
    endpoints.resize(pat.total_spots);
    initpos.resize(pat.total_spots);
    Array4<double> initpos_xbuffer_dbg(pat.total_spots);
    gpu_raytrace_original (pat, ct, endpoints, initpos_xbuffer_dbg, initpos,
                           parser, influence);

    // Warp endpoints in CT ---------------------------
    warped_endpoints.resize(pat.total_spots);
    warped_initpos.resize(pat.total_spots);
    std::copy (endpoints.begin(), endpoints.end(), warped_endpoints.begin());
    std::copy (initpos.begin(), initpos.end(), warped_initpos.begin());
    warper.apply_to(warped_endpoints, warped_initpos, pat.ct, pat.treatment_planes,
                    pat.angles, pat.spots_per_field, parser.warp_opts);

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
                     const Array4<double>& warped_endpoints,
                     const Array4<double>& warped_initpos,
                     std::vector<double>& energy_shift,
                     Array4<double>& cbct_endpoints,
                     Array4<double>& influence)
{
    // Get endpoints in CBCT --------------------
    Volume_t cbct(parser.cbct_file, Volume_t::Source_type::MHA);
    cbct.ext_to_int_coordinates();
    pat.update_geometry_with_external(cbct);

    cbct_endpoints.resize(pat.total_spots);
    gpu_raytrace_warped (pat, cbct, warped_endpoints,
                         warped_initpos, cbct_endpoints,
                         parser, influence);

    // Copy to output vector
    for (size_t i = 0; i < cbct_endpoints.size(); i++)
        energy_shift.at(i) = cbct_endpoints.at(i).w;
    // Apply energy options to output vector
    apply_energy_options(parser.warp_opts, energy_shift, pat.spots_per_field);
}

void export_adapted(Patient_Parameters_t& pat,
                    const Parser& pars,
                    const std::vector<double>& energy_shift,
                    Array4<double>& warped_initpos,
                    Array4<double>& warped_endpoints,
                    const Warper_t& warper)
{
    // Go to virtual source pos
    treatment_plane_to_virtual_src (warped_initpos, warped_endpoints, pat, warper.vf_ave);

    // From virtual source to iso
    virtual_src_to_iso_pos(warped_initpos, pat.virtualSAD);
    
    // Assign to new spotmap
    for (size_t i = 0; i < pat.nbeams; i++)
    {
        std::vector<double>::const_iterator f = energy_shift.begin();
        if (i != 0)
            f += pat.accu_spots_per_field.at(i-1);
        std::vector<double>::const_iterator l = energy_shift.begin() + pat.accu_spots_per_field.at(i);
        std::vector<double> subset_energies(f, l);

        Array4<double>::const_iterator ff = warped_initpos.begin();
        if (i != 0)
            ff += pat.accu_spots_per_field.at(i-1);
        Array4<double>::const_iterator ll = warped_initpos.begin() + pat.accu_spots_per_field.at(i);
        Array4<double> subset_pat_pos(ff, ll);

        if (!pars.output_shifts.empty())
            export_shifts(subset_energies, subset_pat_pos, pars.output_shifts, i, warper.vf_ave);

        Tramp_t tramp(pat.tramp_files.at(i));
        tramp.shift_energies(subset_energies);
        tramp.set_pos(subset_pat_pos);
        tramp.to_file(pat.tramp_files.at(i), pars.out_dir);
    }

}

void export_shifts(const std::vector<double>& e,
                   const Array4<double>& p,
                   const std::string& file,
                   const short& beamid,
                   const Vector3_t<double>& isocenter_shift)
{
    std::ofstream ofs;
    if (beamid != 0)
        ofs.open (file, std::ios::out | std::ios::app);
    else
        ofs.open (file, std::ios::out);

    utils::check_fs(ofs, file, "to write adaptation shifts.");

    if (beamid == 0)
    {
        std::cout << "Writting adaptation shifts to " << file << std::endl;
        std::string txt = "Isocenter shift: " + std::to_string(isocenter_shift.x) +
                          ", " + std::to_string(isocenter_shift.y) +
                          ", " + std::to_string(isocenter_shift.z) + " cm";
        std::cout << txt << std::endl;
        ofs << "# " << txt << '\n';
        ofs << "# e x y z d beamid spotid\n";
    }

    for (size_t spotid = 0; spotid < e.size(); spotid++)
    {
        double vx = p.at(spotid).x;
        double vy = p.at(spotid).y;
        double z = p.at(spotid).z;
        ofs << e.at(spotid) << " " << vx << " " << vy << " " << z << " ";
        ofs << std::sqrt(vx*vx + vy*vy) << " " << beamid << " " << spotid << "\n";
    }
}

void generate_report(const std::string& report,
                     const std::string& output_vf,
                     const std::string& output_shifts,
                     const std::vector<std::string>& tramp_files)
{
    std::cout << "Generating report ..." << std::endl;
    std::string interp = "python3 ";
    std::string code   = std::string(INSTALLATION_PATH) + "/src/extra/create_adapt_report.py ";
    std::string vf     = "--vf " + output_vf + " ";
    std::string shifts = "--shifts " + output_shifts + " ";
    std::string tramps = "--tramps ";
    for (size_t i = 0; i < tramp_files.size(); i++)
        tramps += tramp_files.at(i) + " ";
    std::string outdir = "--outdir " + output_vf.substr(0, output_vf.find_last_of('/')) + " ";
    std::string outfile = "--outfile " + report;
    std::string command = interp + code + vf + shifts + tramps + outdir + outfile;
    std::cout << "Running: " << command << std::endl;
    int res = system(command.c_str());

    if (res)
    {
        std::cerr << "ERROR! Reporting failed with error code: " << res << std::endl;
        exit(EXIT_FAILURE);
    }
}
