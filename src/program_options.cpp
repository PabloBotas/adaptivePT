#include "program_options.hpp"

#include "special_types.hpp"
#include <vector>

void apply_energy_options(Warp_opts_t options,
                          std::vector<float>& energy_shift,
                          const std::vector<short>& spots_per_field)
{
    if (options == FREE_POS_RIGID_ENERGY ||
        options == RIGID_POS_RIGID_ENERGY ||
        options == RIGID_BEAMS_POS_RIGID_ENERGY)
        apply_rigid_energy(energy_shift);
    else if (options == FREE_POS_RIGID_BEAMS_ENERGY ||
             options == RIGID_POS_RIGID_BEAMS_ENERGY ||
             options == RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY)
        apply_rigid_energy_per_beam(energy_shift, spots_per_field);
}

void apply_rigid_energy (std::vector<float>& energy)
{
    std::cout << "Applying rigid energy shifts!!" << std::endl;
    // Calculate average
    float avg = 0;
    for (size_t i = 0; i < energy.size(); i++)
        avg += energy.at(i);
    avg /= energy.size();

    // Set average
    for (size_t i = 0; i < energy.size(); i++)
        energy.at(i) = avg;
}

void apply_rigid_energy_per_beam (std::vector<float>& energy,
                                  const std::vector<short>& spots_per_field)
{
    std::cout << "Applying rigid energy shifts per field!!" << std::endl;
    // Calculate average per field
    float accu_spots = 0;
    std::vector<float> avgs(spots_per_field.size());
    for (size_t ibeam = 0; ibeam < avgs.size(); ibeam++) {
        // std::cout << "BEAM " << ibeam << std::endl;
        for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++) {
            size_t idx;
            if (ibeam > 0)
                idx = ispot + accu_spots;
            else
                idx = ispot;
            avgs.at(ibeam) += energy.at(idx);
            // std::cout << "Index " << idx << " " << energy.at(idx) << " " << avgs.at(ibeam) << std::endl;
        }
        accu_spots += spots_per_field.at(ibeam);
        avgs.at(ibeam) /= spots_per_field.at(ibeam);
        // std::cout << "Average " << avgs.at(ibeam) << std::endl;
    }

    // Set average per field
    accu_spots = 0;
    for (size_t ibeam = 0; ibeam < avgs.size(); ibeam++) {
        for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++) {
            size_t idx;
            if (ibeam > 0)
                idx = ispot+accu_spots;
            else
                idx = ispot;
            energy.at(idx) = avgs.at(ibeam);
        }
        accu_spots += spots_per_field.at(ibeam);
    }
}


// void apply_rigid_layers (const std::vector<short>& spots_per_field,
//                          const std::vector< std::vector<short> >& energy_layers)
// {
//     // Calculate average per layer
//     std::vector<Array3> avgs(spots_per_field.size());
//     for (size_t ibeam = 0; ibeam < energy_layers.size(); ibeam++)
//     {
//         float layer_avg = 0;
//         for (size_t ispot = 0; ispot < energy_layers.at(ibeam); ispot++)
//         {
//             size_t idx = 0;
//             if (ibeam > 0)
//                 idx += energy_layers.at(ibeam-1)
//             layer_avg
//         }
//     }
// }

