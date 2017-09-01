#include <iostream>
#include <iomanip>
#include <string>

#define TO_STRING2(X) #X
#define TO_STRING(X) TO_STRING2(X)
#define INSTALLATION_PATH TO_STRING(BIN_PATH)

#include "command_line_parser.hpp"
#include "initialize_rays.cuh"
#include "gpu_main.cuh"
#include "gpu_source_positioning.cuh"
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

void adjust_pos_in_ct(const Patient_Parameters_t& pat,
                      const Parser& parser,
                      Array4<double>& ct_vf_endpoints,
                      Array4<double>& ct_vf_init_pat_pos,
                      Warper_t& warper);

void adjust_energy_in_cbct(Patient_Parameters_t& pat,
                           const Parser& parser,
                           const Array4<double>& ct_vf_endpoints,
                           const Array4<double>& ct_vf_init_pat_pos,
                           std::vector<double>& energy_shift);

void correct_cold_spots(Patient_Parameters_t& pat,
                        const Parser& parser,
                        const Array4<double>& ct_vf_endpoints,
                        const Array4<double>& ct_vf_init_pat_pos);

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

    // Read input patient
    Patient_Parameters_t pat(parser.patient);
    pat.print();
    pat.ct_to_int_coordinates();
    pat.set_treatment_planes();

    // Start device
    cudaEvent_t start;
    initialize_device(start);

    Array4<double> ct_endpoints(pat.total_spots);
    Array4<double> ct_init_pat_pos(pat.total_spots);
    Warper_t warper(parser.vf_file, parser.output_vf);
    adjust_pos_in_ct (pat, parser, ct_endpoints, ct_init_pat_pos, warper);
    std::vector<double> energy_shift(pat.total_spots);

    adjust_energy_in_cbct (pat, parser, ct_endpoints, ct_init_pat_pos, energy_shift);

    correct_cold_spots (pat, parser.dose_presc, parser.dose_file, ct_endpoints, ct_init_pat_pos);

    export_adapted (pat, parser, energy_shift, ct_init_pat_pos, ct_endpoints, warper);

    if (!parser.report.empty())
        generate_report(parser.report, parser.output_vf, parser.output_shifts, pat.tramp_files);

    // Stop device
    stop_device(start);

    exit(EXIT_SUCCESS);
}

void correct_cold_spots(Patient_Parameters_t& pat,
                        const Parser& parser,
                        const Array4<double>& ct_endpoints,
                        const Array4<double>& ct_init_pat_pos,
                        std::vector<double>& weights)
{
    // Get endpoints in CBCT --------------------
    Volume_t cbct(parser.cbct_file, Volume_t::Source_type::MHA);
    cbct.ext_to_int_coordinates();
    pat.update_geometry_with_external(cbct);

    // Read dose (internal coords) --------------
    Volume_t dose(parser.dose_file, Volume_t::Source_type::DOSE);
    double const MeVg2Gy = 1.6022e-10;
    double const RBE     = 1.1;
    double const flux    = 1.0e9 / parser.spot_factor;
    double scale = RBE*flux*MeVg2Gy;
    dose.normalize(parser.dose_presc/scale);

    gpu_fill_cold_spots (pat, cbct, dose, ct_endpoints, ct_init_pat_pos, weights);
}

void adjust_pos_in_ct(const Patient_Parameters_t& pat,
                      const Parser& parser,
                      Array4<double>& ct_endpoints,
                      Array4<double>& ct_init_pat_pos,
                      Warper_t& warper)
{
    // Read CT and launch rays
    Volume_t ct(pat.planning_ct_file, Volume_t::Source_type::CTVOLUME);
    // The CT volume lacks dimensions information
    ct.setDims(pat.ct);

    // Get endpoints in CT ----------------------------
    ct_endpoints.resize(pat.total_spots);
    ct_init_pat_pos.resize(pat.total_spots);
    Array4<double> ct_init_pos(pat.total_spots);
    gpu_raytrace_original (pat, ct, ct_endpoints, ct_init_pos, ct_init_pat_pos,
                           parser.output_ct_traces);

    // Warp endpoints in CT ---------------------------
    warper.apply_to(ct_endpoints, ct_init_pat_pos, pat.ct, pat.treatment_planes,
                    pat.angles, pat.spots_per_field, parser.warp_opts);

    // Print results
    std::cout << "Warped patient positions and wepl:" << std::endl;
    size_t iters = pat.total_spots < 5 ? pat.total_spots : 5;
    for (size_t i = 0; i < iters; i++)
        ct_init_pat_pos.at(i).print();
}

void adjust_energy_in_cbct(Patient_Parameters_t& pat,
                           const Parser& parser,
                           const Array4<double>& ct_vf_endpoints,
                           const Array4<double>& ct_vf_init_pat_pos,
                           std::vector<double>& energy_shift)
{
    // Get endpoints in CBCT --------------------
    Volume_t cbct(parser.cbct_file, Volume_t::Source_type::MHA);
    cbct.ext_to_int_coordinates();
    pat.update_geometry_with_external(cbct);

    Array4<double> cbct_endpoints(pat.total_spots);
    gpu_raytrace_warped (pat, cbct, ct_vf_endpoints,
                         ct_vf_init_pat_pos, cbct_endpoints,
                         parser.output_cbct_traces);

    // Copy to output vector
    for (size_t i = 0; i < cbct_endpoints.size(); i++)
        energy_shift.at(i) = cbct_endpoints.at(i).w;
    // Apply energy options to output vector
    apply_energy_options(parser.warp_opts, energy_shift, pat.spots_per_field);
}

void export_adapted(Patient_Parameters_t& pat,
                    const Parser& pars,
                    const std::vector<double>& energy_shift,
                    Array4<double>& ct_initpos,
                    Array4<double>& ct_endpoints,
                    const Warper_t& warper)
{
    // Go to virtual source pos
    treatment_plane_to_virtual_src (ct_initpos, ct_endpoints, pat, warper.vf_ave);

    // From virtual source to iso
    virtual_src_to_iso_pos(ct_initpos, pat.virtualSAD);
    
    // Assign to new spotmap
    for (size_t i = 0; i < pat.nbeams; i++)
    {
        std::vector<double>::const_iterator f = energy_shift.begin();
        if (i != 0)
            f += pat.accu_spots_per_field.at(i-1);
        std::vector<double>::const_iterator l = energy_shift.begin() + pat.accu_spots_per_field.at(i);
        std::vector<double> subset_energies(f, l);

        Array4<double>::const_iterator ff = ct_initpos.begin();
        if (i != 0)
            ff += pat.accu_spots_per_field.at(i-1);
        Array4<double>::const_iterator ll = ct_initpos.begin() + pat.accu_spots_per_field.at(i);
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
        std::cout << "Isocenter shift: " << isocenter_shift.x << ", " << isocenter_shift.y << ", " << isocenter_shift.z << std::endl;
        ofs << "# Isocenter shift: " << isocenter_shift.x << " " << isocenter_shift.y << " " << isocenter_shift.z << '\n';
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
