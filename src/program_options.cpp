#include "program_options.hpp"

#include "special_types.hpp"
#include <vector>

void apply_energy_options(Warp_opts_t options,
                          std::vector<double>& energy_shift,
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

void apply_rigid_energy (std::vector<double>& energy)
{
    std::cout << "Applying rigid energy shifts!!" << std::endl;
    // Calculate average
    double avg = 0;
    for (size_t i = 0; i < energy.size(); i++)
        avg += energy.at(i);
    avg /= energy.size();

    // Set average
    for (size_t i = 0; i < energy.size(); i++)
        energy.at(i) = avg;
}

void apply_rigid_energy_per_beam (std::vector<double>& energy,
                                  const std::vector<short>& spots_per_field)
{
    std::cout << "Applying rigid energy shifts per field!!" << std::endl;
    // Calculate average per field
    double accu_spots = 0;
    std::vector<double> avgs(spots_per_field.size());
    for (size_t ibeam = 0; ibeam < avgs.size(); ibeam++)
    {
        // std::cout << "BEAM " << ibeam << std::endl;
        for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++)
        {
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
    for (size_t ibeam = 0; ibeam < avgs.size(); ibeam++)
    {
        for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++)
        {
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

void apply_position_options (Warp_opts_t options,
                             Array3<double>& vf,
                             const std::vector<short>& spots_per_field)
{
    if (options == RIGID_POS_FREE_ENERGY ||
        options == RIGID_POS_RIGID_ENERGY ||
        options == RIGID_POS_RIGID_BEAMS_ENERGY)
        apply_rigid_positions(vf);
    else if (options == RIGID_BEAMS_POS_FREE_ENERGY ||
             options == RIGID_BEAMS_POS_RIGID_ENERGY ||
             options == RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY)
        apply_rigid_positions_per_beam(vf, spots_per_field);
    
}

void apply_rigid_positions (Array3<double>& vf)
{
    std::cout << "Applying rigid positions!!" << std::endl;
    // Calculate average
    Vector3_t<double> avg;
    for (size_t i = 0; i < vf.size(); i++)
        avg += vf.at(i);
    avg /= vf.size();

    // Set average
    for (size_t i = 0; i < vf.size(); i++)
        vf.at(i) = avg;
}

void apply_rigid_positions_per_beam (Array3<double>& vf,
                                     const std::vector<short>& spots_per_field)
{
    std::cout << "Applying rigid positions per field!!" << std::endl;
    // Calculate average per field
    double accu_spots = 0;
    Array3<double> avgs(spots_per_field.size());
    for (size_t ibeam = 0; ibeam < avgs.size(); ibeam++)
    {
        // std::cout << "BEAM " << ibeam << std::endl;
        for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++)
        {
            size_t idx;
            if (ibeam > 0)
                idx = ispot + accu_spots;
            else
                idx = ispot;
            avgs.at(ibeam) += vf.at(idx);
            // std::cout << "Index " << idx << " " << vf.at(idx).x << " " << avgs.at(ibeam).x << std::endl;
        }
        accu_spots += spots_per_field.at(ibeam);
        avgs.at(ibeam) /= spots_per_field.at(ibeam);
        // std::cout << "Average " << avgs.at(ibeam).x << std::endl;
    }

    // Set average per field
    accu_spots = 0;
    for (size_t ibeam = 0; ibeam < avgs.size(); ibeam++)
    {
        for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++)
        {
            size_t idx;
            if (ibeam > 0)
                idx = ispot+accu_spots;
            else
                idx = ispot;
            vf.at(idx) = avgs.at(ibeam);
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
//         double layer_avg = 0;
//         for (size_t ispot = 0; ispot < energy_layers.at(ibeam); ispot++)
//         {
//             size_t idx = 0;
//             if (ibeam > 0)
//                 idx += energy_layers.at(ibeam-1)
//             layer_avg
//         }
//     }
// }

