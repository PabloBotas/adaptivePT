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

void deal_with_ct(const Patient_Parameters_t& pat,
                  const Parser& parser,
                  Array4<double>& ct_vf_endpoints,
                  Array4<double>& ct_vf_init_pat_pos);

void deal_with_cbct(Patient_Parameters_t& pat,
                    const Parser& parser,
                    const Array4<double>& ct_vf_endpoints,
                    const Array4<double>& ct_vf_init_pat_pos,
                    std::vector<double>& energy_shift);

void export_adapted(Patient_Parameters_t& pat,
                    const Parser& pars,
                    const std::vector<double>& energy_shift,
                    Array4<double>& pat_pos,
                    Array4<double>& pat_pos2);

void export_shifts(const std::vector<double>& e,
                   const Array4<double>& p,
                   const std::string& file,
                   const short& beamid);

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
    deal_with_ct (pat, parser, ct_endpoints, ct_init_pat_pos);
    std::vector<double> energy_shift(pat.total_spots);
    if (!parser.no_energy)
        deal_with_cbct (pat, parser, ct_endpoints, ct_init_pat_pos, energy_shift);

    export_adapted (pat, parser, energy_shift, ct_init_pat_pos, ct_endpoints);

    if (!parser.report.empty())
        generate_report(parser.report, parser.output_vf, parser.output_shifts, pat.tramp_files);

    // Stop device
    stop_device(start);

    exit(EXIT_SUCCESS);
}

void deal_with_ct(const Patient_Parameters_t& pat,
                  const Parser& parser,
                  Array4<double>& ct_endpoints,
                  Array4<double>& ct_init_pat_pos)
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
    Warper_t warp(parser.vf_file, parser.output_vf);
    warp.apply_to(ct_endpoints, ct_init_pat_pos, pat.ct, pat.treatment_planes,
                  pat.spots_per_field, pat.angles);

    // Print results
    std::cout << "Warped patient positions and wepl:" << std::endl;
    size_t iters = pat.total_spots < 5 ? pat.total_spots : 5;
    for (size_t i = 0; i < iters; i++)
        ct_init_pat_pos.at(i).print();
}

void deal_with_cbct(Patient_Parameters_t& pat,
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

    // Print results ----------------------------
    for (size_t i = 0; i < cbct_endpoints.size(); i++)
        energy_shift.at(i) = cbct_endpoints.at(i).w;
}

void export_adapted(Patient_Parameters_t& pat,
                    const Parser& pars,
                    const std::vector<double>& energy_shift,
                    Array4<double>& pat_pos,
                    Array4<double>& pat_pos2)
{
    // Go to virtual source pos
    treatment_plane_to_virtual_src (pat_pos, pat_pos2, pat);

    // From virtual source to iso
    virtual_src_to_iso_pos(pat_pos, pat.virtualSAD);
    
    // Assign to new spotmap
    for (size_t i = 0; i < pat.nbeams; i++)
    {
        std::vector<double>::const_iterator f = energy_shift.begin();
        if (i != 0)
            f += pat.accu_spots_per_field.at(i-1);
        std::vector<double>::const_iterator l = energy_shift.begin() + pat.accu_spots_per_field.at(i);
        std::vector<double> subset_energies(f, l);

        Array4<double>::const_iterator ff = pat_pos.begin();
        if (i != 0)
            ff += pat.accu_spots_per_field.at(i-1);
        Array4<double>::const_iterator ll = pat_pos.begin() + pat.accu_spots_per_field.at(i);
        Array4<double> subset_pat_pos(ff, ll);

        Tramp_t tramp(pat.tramp_files.at(i));

        if (!pars.output_shifts.empty())
            export_shifts(subset_energies, subset_pat_pos, pars.output_shifts, i);
        if (!pars.no_energy)
            tramp.shift_energies(subset_energies);
        tramp.set_pos(subset_pat_pos);
        tramp.to_file(pat.tramp_files.at(i), pars.out_dir);
    }
}

void export_shifts(const std::vector<double>& e,
                   const Array4<double>& p,
                   const std::string& file,
                   const short& beamid)
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
        ofs << "e x y z d beamid spotid\n";
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
