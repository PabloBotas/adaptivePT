#include <iostream>
#include <string>

#include "command_line_parser.hpp"
#include "patient_parameters.hpp"
#include "volume.hpp"
#include "gpu_main.cuh"

// TODO DONE Read CT
// TODO Get endpoints on CT
// TODO Apply VF to endpoints -> updated endpoints
// TODO Apply VF to tramp (2D) <- geometrical adjustment
// TODO DONE Read CBCT
// TODO Get endpoints on CBCT of updated tramp files
// TODO Get energy shift between CBCT endpoints and update CT endpoints
// TODO Correct tramp again (E) <- energy adjustment

int main(int argc, char** argv)
{
    Parser parser(argc, argv);
    parser.print_parameters();

    // Read input patient
    Patient_Parameters_t patient_data(parser.patient);
    patient_data.add_results_directory(parser.out_dir);
    patient_data.print();
    patient_data.adjust_to_internal_coordinates();

    // Start device
    cudaEvent_t start, stop;
    initialize_device(start, stop);

    // Read CT and launch rays
    Patient_Volume_t ct(patient_data.planning_ct_file, Patient_Volume_t::Source_type::CTVOLUME);
    // The CT volume lacks dimensions information

    ct.setDimsFromPatient(patient_data.ct);
    // Get endpoints
    std::vector<float4> ct_endpoints = gpu_get_beam_endpoints(patient_data, ct);

    // Print results
    std::cout << "SpotID \t WEPL \t X \t Y \t Z \t WEPL" << std::endl;
    for (size_t i = 0; i < 10; i++)
    {
        std::cout << i << "\t" << ct_endpoints.at(i).w;
        std::cout << "\t" << ct_endpoints.at(i).x;
        std::cout << "\t" << ct_endpoints.at(i).y;
        std::cout << "\t" << ct_endpoints.at(i).z;
        std::cout << "\t" << ct_endpoints.at(i).w << std::endl;
    }

    // Stop device
    stop_device(start, stop);

//    // Read CBCT and launch rays
//    Patient_Volume_t cbct(parser.cbct_file, Patient_Volume_t::Source_type::MHA);
//    patient_data.update_geometry_offsets(cbct);
//    std::vector<float4> cbct_endpoints = gpu_get_beam_endpoints(patient_data, cbct);

    // Get geometric coordinates of such WEPLS


    // Apply vector field to those coordinates (same applied to contours)
    // Get WEPLS of the coordinates in CBCT
    // Adapt initial energies to new geometry

    // plan_in.print(0);
    // plan_in.to_file("test.tramp");

    // Finalize the entire computation

    exit(EXIT_SUCCESS);
}


