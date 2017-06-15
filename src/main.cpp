#include <iostream>
#include <string>

#include "command_line_parser.hpp"
#include "gpu_main.cuh"
#include "gpu_source_positioning.cuh"
#include "patient_parameters.hpp"
#include "utils.hpp"
#include "volume.hpp"
#include "warper.hpp"


// TODO DONE Read CT
// TODO DONE Get endpoints on CT
// TODO DONE Apply VF to endpoints -> updated endpoints
// TODO DONE Read CBCT
// TODO DONE Projection onto treatment plane
// TODO Backtrace vf_endpoints in CBCT until initial wepl
// TODO Get distances to initial endpoint positions
// TODO Get WEPL and energy distance

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

    // Read CT and launch rays
    Patient_Volume_t ct(patient_data.planning_ct_file, Patient_Volume_t::Source_type::CTVOLUME);
    // The CT volume lacks dimensions information
    ct.setDims(patient_data.ct);

    // Get endpoints in CT ----------------------------
    std::vector< Vector4_t<float> > ct_endpoints(patient_data.total_spots);
    std::vector< Vector4_t<float> > ct_init_pat_pos(patient_data.total_spots);
    std::vector< Vector4_t<float> > ct_init_pos(patient_data.total_spots);
    gpu_raytrace_original (patient_data, ct, ct_endpoints, ct_init_pos, ct_init_pat_pos,
                           "output_volume.raw");
    std::vector< Vector4_t<float> > treatment_plane = get_treatment_planes(patient_data.angles);
    // Print results
    size_t iters = ct_endpoints.size() < 5 ? ct_endpoints.size() : 5;

    // Warp endpoints in CT ---------------------------
    std::vector< Vector4_t<float> > ct_vf_endpoints = ct_endpoints;
    std::vector< Vector4_t<float> > ct_vf_init_pat_pos = ct_init_pat_pos;
    warp_data (ct_vf_endpoints, ct_vf_init_pat_pos, parser.vf_file,
               patient_data.ct, treatment_plane);
    
    // Print results
    std::cout << "\nPatient positions and wepl:" << std::endl;
    for (size_t i = 0; i < iters; i++)
        ct_init_pat_pos.at(i).print();
    std::cout << "Warped patient positions and wepl:" << std::endl;
    for (size_t i = 0; i < iters; i++)
        ct_vf_init_pat_pos.at(i).print();

    // Get endpoints in CBCT --------------------
    std::vector< Vector4_t<float> > cbct_endpoints(patient_data.total_spots, 5);
    gpu_raytrace_warped (patient_data, ct, ct_vf_endpoints,
                         ct_vf_init_pat_pos, cbct_endpoints,
                         "output_volume_cbct.raw");

    // Print results
    std::cout << "X \t Y \t Z \t Energy adapt" << std::endl;
    for (size_t i = 0; i < iters; i++)
        cbct_endpoints.at(i).print();
    
    // std::cout << "X \t Y \t Z \t E" << std::endl;

    // Stop device
    stop_device(start);

//    // Read CBCT and launch rays
//    Patient_Volume_t cbct(parser.cbct_file, Patient_Volume_t::Source_type::MHA);
//    patient_data.update_geometry_offsets(cbct);
//    std::vector<float4> cbct_endpoints = gpu_get_beam_endpoints(patient_data, cbct);

    // Finalize the entire computation

    exit(EXIT_SUCCESS);
}


