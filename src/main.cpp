#include <iostream>
#include <string>

#include "command_line_parser.hpp"
#include "tramp.hpp"

int main(int argc, char** argv)
{
    Parser parser(argc, argv);
    parser.print_parameters();

    // Read input tramp
    Tramp_t tramp_in(parser.tramp_in_file);

    // Calculate WEPL from energies
    std::vector<double> wepls = tramp_in.getWEPLs();
    for (size_t i = 0; i < tramp_in.nspots; i++)
    {
        double energy = tramp_in.spots.at(i).e;
        std::cout << energy << " " << wepls.at(i) << std::endl;
    }

    // Get geometric coordinates of such WEPLS
    // Apply vector field to those coordinates (same applied too contours)
    // Get WEPLS of the coordinates in CBCT
    // Adapt initial energies to new geometry

    // tramp_in.print(0);
    tramp_in.to_file("test.tramp");
    
    exit(EXIT_SUCCESS);
}






