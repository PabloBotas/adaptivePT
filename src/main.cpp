#include <iostream>
#include <string>

#include "command_line_parser.hpp"
#include "gpu_main.cuh"
#include "patient_parameters.hpp"
#include "utils.hpp"
#include "volume.hpp"
#include "warper.hpp"


// TODO DONE Read CT
// TODO DONE Get endpoints on CT
// TODO DONE Apply VF to endpoints -> updated endpoints
// TODO DONE Read CBCT
// TODO Projection onto treatment plane is done per spot,
//      which is not correct because they are not parallel to each other: FIX
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
    std::vector< Vector4_t<float> > ct_init_pos(patient_data.total_spots);
    std::vector< Vector2_t<short> > ct_metadata(patient_data.total_spots);
    gpu_raytrace_plan(patient_data, ct, ct_endpoints, ct_init_pos, ct_metadata);
    // Print results
    size_t iters = ct_endpoints.size() < 5 ? ct_endpoints.size() : 5;
    std::cout << "X \t Y \t Z \t WEPL" << std::endl;
    for (size_t i = 0; i < iters; i++)
        ct_endpoints.at(i).print();

    // Warp endpoints in CT ---------------------------
    std::vector< Vector4_t<float> > ct_vf_endpoints = ct_endpoints;
    std::vector< Vector4_t<float> > ct_vf_init_pos = ct_init_pos;
    warp_data(ct_vf_endpoints, ct_vf_init_pos, parser.vf_file, patient_data.ct);
    
    // Print results
    std::cout << "X \t Y \t Z \t WEPL" << std::endl;
    for (size_t i = 0; i < iters; i++)
        ct_vf_endpoints.at(i).print();

    // Backtrace endpoints in CBCT --------------------
    std::vector< Vector4_t<float> > adapted_sources;
    gpu_backtrace_endpoints(patient_data, ct,
                            ct_vf_endpoints, ct_init_pos, ct_metadata,
                            adapted_sources);

    // Stop device
    stop_device(start);

//    // Read CBCT and launch rays
//    Patient_Volume_t cbct(parser.cbct_file, Patient_Volume_t::Source_type::MHA);
//    patient_data.update_geometry_offsets(cbct);
//    std::vector<float4> cbct_endpoints = gpu_get_beam_endpoints(patient_data, cbct);

    // Finalize the entire computation

    exit(EXIT_SUCCESS);
}


