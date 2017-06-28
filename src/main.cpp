#include <iostream>
#include <string>

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
// TODO Convert data back to tramp file
// TODO Fix mha output's offset to visualize in Slicer
// TODO Verify that VF probing is correctly done. How?
// Done Verify CBCT offsets are correctly updated

void deal_with_ct(Patient_Parameters_t& pat,
                  const Parser& parser,
                  Array4<float>& ct_vf_endpoints,
                  Array4<float>& ct_vf_init_pat_pos);

void deal_with_cbct(Patient_Parameters_t& pat,
                    const Parser& parser,
                    const Array4<float>& ct_vf_endpoints,
                    const Array4<float>& ct_vf_init_pat_pos,
                    std::vector<float>& energy_shift);

void export_adapted(Patient_Parameters_t& pat,
                    const std::string& out_dir,
                    const std::vector<float>& energy_shift,
                    Array4<float>& pat_pos,
                    Array4<float>& pat_pos2);

int main(int argc, char** argv)
{
    Parser parser(argc, argv);
    parser.print_parameters();

    // Read input patient
    Patient_Parameters_t pat(parser.patient);
    pat.add_results_directory(parser.out_dir);
    pat.print();
    pat.ext_to_int_coordinates();
    pat.set_treatment_planes();

    // Start device
    cudaEvent_t start;
    initialize_device(start);

    Array4<float> ct_endpoints(pat.total_spots);
    Array4<float> ct_init_pat_pos(pat.total_spots);
    deal_with_ct (pat, parser, ct_endpoints, ct_init_pat_pos);
    std::vector<float> energy_shift(pat.total_spots);
    deal_with_cbct (pat, parser, ct_endpoints, ct_init_pat_pos, energy_shift);

    export_adapted (pat, parser.out_dir, energy_shift, ct_init_pat_pos, ct_endpoints);

    // Stop device
    stop_device(start);

    exit(EXIT_SUCCESS);
}

void deal_with_ct(Patient_Parameters_t& pat,
                  const Parser& parser,
                  Array4<float>& ct_endpoints,
                  Array4<float>& ct_init_pat_pos)
{
    // Read CT and launch rays
    Volume_t ct(pat.planning_ct_file, Volume_t::Source_type::CTVOLUME);
    // The CT volume lacks dimensions information
    ct.setDims(pat.ct);

    // Get endpoints in CT ----------------------------
    ct_endpoints.resize(pat.total_spots);
    ct_init_pat_pos.resize(pat.total_spots);
    Array4<float> ct_init_pos(pat.total_spots);
    gpu_raytrace_original (pat, ct, ct_endpoints, ct_init_pos, ct_init_pat_pos,
                           "output_volume.raw");
    // Print results
    size_t iters = pat.total_spots < 5 ? pat.total_spots : 5;

    // Warp endpoints in CT ---------------------------
    warp_data (ct_endpoints, ct_init_pat_pos, parser.vf_file,
               pat.ct, pat.treatment_planes.dir, pat.spots_per_field);
    
    // Print results
    std::cout << "Warped patient positions and wepl:" << std::endl;
    for (size_t i = 0; i < iters; i++)
        ct_init_pat_pos.at(i).print();
}

void deal_with_cbct(Patient_Parameters_t& pat,
                    const Parser& parser,
                    const Array4<float>& ct_vf_endpoints,
                    const Array4<float>& ct_vf_init_pat_pos,
                    std::vector<float>& energy_shift)
{
    // Get endpoints in CBCT --------------------
    Volume_t cbct(parser.cbct_file, Volume_t::Source_type::MHA);
    cbct.ext_to_int_coordinates();
    pat.update_geometry_with_external(cbct);

    Array4<float> cbct_endpoints(pat.total_spots);
    gpu_raytrace_warped (pat, cbct, ct_vf_endpoints,
                         ct_vf_init_pat_pos, cbct_endpoints,
                         "output_volume_cbct.raw");

    // Print results ----------------------------
    for (size_t i = 0; i < cbct_endpoints.size(); i++)
        energy_shift.at(i) = cbct_endpoints.at(i).w;
}

void export_adapted(Patient_Parameters_t& pat,
                    const std::string& out_dir,
                    const std::vector<float>& energy_shift,
                    Array4<float>& pat_pos,
                    Array4<float>& pat_pos2)
{
    // Go to virtual source pos
    treatment_plane_to_virtual_src (pat_pos, pat_pos2, pat);

    // From virtual source to iso
    virtual_src_to_iso_pos(pat_pos, pat.virtualSAD);
    
    // Assign to new spotmap
    for (size_t i = 0; i < pat.nbeams; i++)
    {
        Tramp_t tramp(pat.tramp_files.at(i));
        tramp.shift_energies(energy_shift);
        tramp.set_pos(pat_pos);
        tramp.to_file(pat.tramp_files.at(i), out_dir);
    }
}
