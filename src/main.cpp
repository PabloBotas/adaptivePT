#include <iostream>
#include <string>

#include "command_line_parser.hpp"
#include "patient_parameters.hpp"
#include "volume.hpp"

int main(int argc, char** argv)
{
    Parser parser(argc, argv);
    parser.print_parameters();

    // Read input patient
    Patient_Parameters_t patient_data(parser.patient);
    patient_data.print();

    // Read CBCT
    Patient_Volume_t cbct(parser.cbct_file, Patient_Volume_t::Source_type::MHA_CALIBRATION);
    std::cout << cbct.file << std::endl;
    std::cout << cbct.n.x << std::endl;
    std::cout << cbct.n.y << std::endl;
    std::cout << cbct.n.z << std::endl;
    std::cout << cbct.d.x << std::endl;
    std::cout << cbct.d.y << std::endl;
    std::cout << cbct.d.z << std::endl;
    std::cout << cbct.nElements << std::endl;
    float m=0;
    for (size_t i = 0; i < cbct.nElements; i++)
    {
        m = cbct.hu.at(i) > m ? cbct.hu.at(i) : m;
    }
    std::cout << "Maximum: " << m << std::endl;


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






