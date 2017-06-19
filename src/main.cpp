#include <iostream>
#include <string>

#include "command_line_parser.hpp"
#include "gpu_main.cuh"
#include "gpu_source_positioning.cuh"
#include "patient_parameters.hpp"
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
// TODO Verify CBCT offsets are correctly updated

void deal_with_ct(Patient_Parameters_t& patient_data,
                  const Parser& parser,
                  Array4<float>& ct_vf_endpoints,
                  Array4<float>& ct_vf_init_pat_pos);

void deal_with_cbct(Patient_Parameters_t& patient_data,
                    const Parser& parser,
                    Array4<float>& ct_vf_endpoints,
                    Array4<float>& ct_vf_init_pat_pos);

int main(int argc, char** argv)
{
    Parser parser(argc, argv);
    parser.print_parameters();

    // Read input patient
    Patient_Parameters_t patient_data(parser.patient);
    patient_data.add_results_directory(parser.out_dir);
    patient_data.print();
    patient_data.ext_to_int_coordinates();

    // Start device
    cudaEvent_t start;
    initialize_device(start);

    Array4<float> ct_endpoints(patient_data.total_spots);
    Array4<float> ct_init_pat_pos(patient_data.total_spots);
    deal_with_ct (patient_data, parser, ct_endpoints, ct_init_pat_pos);
    deal_with_cbct (patient_data, parser, ct_endpoints, ct_init_pat_pos);

    // Stop device
    stop_device(start);

    exit(EXIT_SUCCESS);
}

void deal_with_ct(Patient_Parameters_t& patient_data,
                  const Parser& parser,
                  Array4<float>& ct_endpoints,
                  Array4<float>& ct_init_pat_pos)
{
    // Read CT and launch rays
    Patient_Volume_t ct(patient_data.planning_ct_file, Patient_Volume_t::Source_type::CTVOLUME);
    // The CT volume lacks dimensions information
    ct.setDims(patient_data.ct);

    // Get endpoints in CT ----------------------------
    ct_endpoints.resize(patient_data.total_spots);
    ct_init_pat_pos.resize(patient_data.total_spots);
    Array4<float> ct_init_pos(patient_data.total_spots);
    gpu_raytrace_original (patient_data, ct, ct_endpoints, ct_init_pos, ct_init_pat_pos,
                           "output_volume.raw");
    Array4<float> treatment_plane = get_treatment_planes(patient_data.angles);
    // Print results
    size_t iters = patient_data.total_spots < 5 ? patient_data.total_spots : 5;

    // Warp endpoints in CT ---------------------------
    warp_data (ct_endpoints, ct_init_pat_pos, parser.vf_file,
               patient_data.ct, treatment_plane);
    
    // Print results
    std::cout << "Warped patient positions and wepl:" << std::endl;
    for (size_t i = 0; i < iters; i++)
        ct_init_pat_pos.at(i).print();
}

void deal_with_cbct(Patient_Parameters_t& patient_data,
                    const Parser& parser,
                    Array4<float>& ct_vf_endpoints,
                    Array4<float>& ct_vf_init_pat_pos)
{
    // Get endpoints in CBCT --------------------
    Patient_Volume_t cbct(parser.cbct_file, Patient_Volume_t::Source_type::MHA);
    patient_data.update_geometry_offsets(cbct);

    Array4<float> cbct_endpoints(patient_data.total_spots);
    gpu_raytrace_warped (patient_data, cbct, ct_vf_endpoints,
                         ct_vf_init_pat_pos, cbct_endpoints,
                         "output_volume_cbct.raw");

    // Print results
    size_t iters = patient_data.total_spots < 5 ? patient_data.total_spots : 5;
    std::cout << "X \t Y \t Z \t Energy adapt" << std::endl;
    for (size_t i = 0; i < iters; i++)
        cbct_endpoints.at(i).print();
}
