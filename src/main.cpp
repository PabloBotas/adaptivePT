#include <iostream>
#include <string>

#include "command_line_parser.hpp"
#include "patient_parameters.hpp"
#include "volume.hpp"
#include "gpu_main.cuh"

#include "utils.hpp"

// TODO DONE Read CT
// TODO DONE Get endpoints on CT
// TODO DONE Apply VF to endpoints -> updated endpoints
// TODO DONE Read CBCT
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
    cudaEvent_t start, stop;
    initialize_device(start, stop);

    // Read CT and launch rays
    Patient_Volume_t ct(patient_data.planning_ct_file, Patient_Volume_t::Source_type::CTVOLUME);
    // The CT volume lacks dimensions information
    ct.setDims(patient_data.ct);
    // Get endpoints
    std::vector< Vector4_t<float> > ct_endpoints = gpu_raytrace_plan(patient_data, ct);

    // Print results
    size_t iters = ct_endpoints.size() < 10 ? ct_endpoints.size() : 10;
    std::cout << "X \t Y \t Z \t WEPL" << std::endl;
    for (size_t i = 0; i < iters; i++)
        ct_endpoints.at(i).print();

    utils::flip_positions_X(ct_endpoints, patient_data.ct);
    std::vector< Vector4_t<float> > ct_vf_endpoints = ct_endpoints;
    utils::run_plastimatch_probe(ct_vf_endpoints, parser.vf_file);

    std::cout << "X \t Y \t Z \t WEPL" << std::endl;
    for (size_t i = 0; i < iters; i++)
        ct_endpoints.at(i).print();

    // Stop device
    stop_device(start, stop);

//    // Read CBCT and launch rays
//    Patient_Volume_t cbct(parser.cbct_file, Patient_Volume_t::Source_type::MHA);
//    patient_data.update_geometry_offsets(cbct);
//    std::vector<float4> cbct_endpoints = gpu_get_beam_endpoints(patient_data, cbct);

    // Finalize the entire computation

    exit(EXIT_SUCCESS);
}


