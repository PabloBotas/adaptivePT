#include <iostream>
#include <string>

#include "command_line_parser.hpp"
#include "plan_parameters.hpp"
#include "mha_volumes.hpp"

int main(int argc, char** argv)
{
    Parser parser(argc, argv);
    parser.print_parameters();

    // Read input patient
    Plan_Parameters_t patient_data(parser.patient);
    patient_data.print();

    // Read CBCT
    Mha_t cbct(parser.cbct_file);


    // Calculate WEPL from energies
    // std::vector<double> wepls = plan_in.getWEPLs();

    // Get geometric coordinates of such WEPLS


    // Apply vector field to those coordinates (same applied too contours)
    // Get WEPLS of the coordinates in CBCT
    // Adapt initial energies to new geometry

    // plan_in.print(0);
    // plan_in.to_file("test.tramp");
    
    exit(EXIT_SUCCESS);
}






