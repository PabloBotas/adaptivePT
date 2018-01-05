#include "program_options.hpp"

#include "special_types.hpp"
#include "utils.hpp"
#include <algorithm>
#include <vector>

void calculate_range_shifter (const RShifter_steps_t& rshifter_steps,
                              const Adapt_constraints_t& constr,
                              std::vector<float>& new_energies,
                              std::vector<RangeShifter_Dims_t>& rs,
                              std::string& machine,
                              const std::vector<float>& isocenter_to_beam_distance,
                              const std::vector<float>& orig_energies,
                              const std::vector<short>& spots_per_field)
{
    if (constr == Adapt_constraints_t::RANGE_SHIFTER ||
        constr == Adapt_constraints_t::ISOCENTER_SHIFT_RANGE_SHIFTER ||
        constr == Adapt_constraints_t::FIELD_ISOCENTER_SHIFT_RANGE_SHIFTER) {
        physical_range_shifter (rshifter_steps, new_energies, rs, isocenter_to_beam_distance, 
                                orig_energies, spots_per_field);
    } else if (constr == Adapt_constraints_t::V_RANGE_SHIFTER ||
               constr == Adapt_constraints_t::ISOCENTER_SHIFT_V_RANGE_SHIFTER ||
               constr == Adapt_constraints_t::FIELD_ISOCENTER_SHIFT_V_RANGE_SHIFTER) {
        virtual_range_shifter (rshifter_steps, new_energies, orig_energies, spots_per_field);
    } else {
        std::vector<float> below_thres;
        if (outside_machine_energies(below_thres, new_energies, machine, spots_per_field)) {
            correct_energy_range (rshifter_steps, new_energies, rs,
                                  isocenter_to_beam_distance, 
                                  spots_per_field, below_thres, machine);
        }
    } 
}


bool outside_machine_energies (std::vector<float>& below_thres,
                               std::vector<float>& energies,
                               const std::string& machine,
                               const std::vector<short>& spots_per_field)
{
    float min_e{}, max_e{};
    if (machine == "topassmallSpots"     ||
        machine == "topas_a5_smallSpots" ||
        machine == "topasmediumspots"    ||
        machine == "topas_a5_mediumspots") {
        min_e = 59.769*1e6;
        max_e = 230.527*1e6;
    } else if (machine == "topasmghr4" ||
               machine == "topasmghr5") {
        min_e = 91.015*1e6;
        max_e = 223.58*1e6;
    }

    uint nbeams = spots_per_field.size();
    below_thres.resize(nbeams, max_e);
    bool if_need_correction = false;
    uint accu_spots{};
    for (size_t ibeam = 0; ibeam < nbeams; ibeam++) {
        for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++) {
            uint idx = ispot + accu_spots;
            if (energies.at(idx) > max_e) {
                std::cerr << "WARNING! Energy of spot @ beam: " << ispot << " @ " << ibeam;
                std::cerr << " = " << energies.at(idx)/1e6 << " MeV. ";
                std::cerr << "Capped to " << machine << " gantry max energy of ";
                std::cerr << max_e/1e6 << " MeV." << std::endl;
                energies.at(idx) = max_e;
            }
            if (energies.at(idx) < min_e) {
                std::cerr << "WARNING! Energy of spot @ beam: " << ispot << " @ " << ibeam;
                std::cerr << " = " << energies.at(idx)/1e6 << " MeV < " << min_e/1e6;
                std::cerr << " MeV (gantry min). We will attempt to increase it and add a range "
                             "shifter." << std::endl;
                below_thres.at(ibeam) = std::min(energies.at(idx), below_thres.at(ibeam));
                if_need_correction = true;
            }
        }
        if (below_thres.at(ibeam) == max_e)
            below_thres.at(ibeam) = 0.f;
        accu_spots += spots_per_field.at(ibeam);
    }

    return if_need_correction;
}



void virtual_range_shifter (const RShifter_steps_t& rshifter_steps,
                            std::vector<float>& new_energies,
                            const std::vector<float>& orig_energies,
                            const std::vector<short>& spots_per_field)
{
    std::cout << "Applying virtual range shifter per field!!" << std::endl;
    // Calculate average per field
    uint nbeams = spots_per_field.size();
    uint accu_spots = 0;
    std::vector<float> range_change_ave(nbeams, 0);

    const double e_r[6] = {-0.04796, 0.008392, 0.0008097, -1.357e-06, 2.099e-09, -1.754e-12};
    const double r_e[6] = {7.99143832, 6.69077816e-01, 1.96456917e-03,
                           1.80100513, 7.32441249e-08, -2.11977495e-01};

    for (size_t i = 0; i < nbeams; i++) {
        // Average energy of all the spots
        float average_energy = 
            std::accumulate(orig_energies.begin(), orig_energies.end(), 0.0)/orig_energies.size();

        // Get average range change of the spots with initial E >= average_energy
        uint considered_spots = 0;
        for (short ispot = 0; ispot < spots_per_field.at(i); ispot++) {
            if (orig_energies.at(ispot + accu_spots) >= average_energy) {
                const float& E1 = new_energies.at(ispot + accu_spots);
                const float& E0 = orig_energies.at(ispot + accu_spots);
                const float R1 = e_r[0] + e_r[1]*E1/1e6 + e_r[2]*std::pow(E1/1e6, 2) +
                                 e_r[3]*std::pow(E1/1e6, 3) + e_r[4]*std::pow(E1/1e6, 4) +
                                 e_r[5]*std::pow(E1/1e6, 5);
                const float R0 = e_r[0] + e_r[1]*E0 + e_r[2]*std::pow(E0, 2) +
                                 e_r[3]*std::pow(E0, 3) + e_r[4]*std::pow(E0, 4) +
                                 e_r[5]*std::pow(E0, 5);
                range_change_ave.at(i) += R1-R0;
                considered_spots++;
            }
        }
        range_change_ave.at(i) /= considered_spots;
        range_change_ave.at(i) = rshifter_steps.select_range_shifter(range_change_ave.at(i));
        accu_spots += spots_per_field.at(i);
    }

    accu_spots = 0;
    for (size_t i = 0; i < nbeams; i++) {
        if (range_change_ave.at(i) != 0) {
            // what is the energy shift per layer so that we get the remaining range shift?
            float E_new;
            for (short ispot = 0; ispot < spots_per_field.at(i); ispot++) {
                uint idx = ispot + accu_spots;
                if (!(ispot != 0 && orig_energies.at(idx) == orig_energies.at(idx-1))) {
                    const float& E0 = orig_energies.at(idx);
                    float R0 = e_r[0] + e_r[1]*E0 + e_r[2]*std::pow(E0, 2) +
                               e_r[3]*std::pow(E0, 3) + e_r[4]*std::pow(E0, 4) +
                               e_r[5]*std::pow(E0, 5);
                    float R = R0 + range_change_ave.at(i);
                    E_new = std::pow(R/r_e[0], 1/r_e[1]) + std::pow(R/r_e[2], 1/r_e[3]) +
                            std::pow(R/r_e[4], 1/r_e[5]);
                }
                new_energies.at(idx) = E_new*1e6;
            }
        } else {
            for (short ispot = 0; ispot < spots_per_field.at(i); ispot++) {
                new_energies.at(ispot + accu_spots) = orig_energies.at(ispot + accu_spots);
            }
        }
        accu_spots += spots_per_field.at(i);
    }
}


void physical_range_shifter (const RShifter_steps_t& rshifter_steps,
                             std::vector<float>& new_energies,
                             std::vector<RangeShifter_Dims_t>& rs,
                             const std::vector<float>& isocenter_to_beam_distance,
                             const std::vector<float>& orig_energies,
                             const std::vector<short>& spots_per_field)
{
    std::cout << "Applying physical range shifter per field!!" << std::endl;
    // Calculate average per field
    uint nbeams = spots_per_field.size();
    uint accu_spots = 0;
    std::vector<float> range_change_ave(nbeams);

    const double e_r[6] = {-0.04796, 0.008392, 0.0008097, -1.357e-06, 2.099e-09, -1.754e-12};
    const double r_e[6] = {7.99143832, 6.69077816e-01, 1.96456917e-03,
                           1.80100513, 7.32441249e-08, -2.11977495e-01};
    for (size_t ibeam = 0; ibeam < nbeams; ibeam++) {
        float average_energy = 0;
        for (uint i = 0; i < orig_energies.size(); ++i)
            average_energy += orig_energies.at(i);
        average_energy /= orig_energies.size();

        uint considered_spots = 0;
        for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++) {
            if (orig_energies.at(ispot + accu_spots) >= average_energy) {
                const float& E1 = new_energies.at(ispot + accu_spots);
                const float& E0 = orig_energies.at(ispot + accu_spots);
                const float R1 = e_r[0] + e_r[1]*E1/1e6 + e_r[2]*std::pow(E1/1e6, 2) +
                                 e_r[3]*std::pow(E1/1e6, 3) + e_r[4]*std::pow(E1/1e6, 4) +
                                 e_r[5]*std::pow(E1/1e6, 5);
                const float R0 = e_r[0] + e_r[1]*E0 + e_r[2]*std::pow(E0, 2) +
                                 e_r[3]*std::pow(E0, 3) + e_r[4]*std::pow(E0, 4) +
                                 e_r[5]*std::pow(E0, 5);
                range_change_ave.at(ibeam) += R1-R0;
                considered_spots++;
            }
        }
        range_change_ave.at(ibeam) /= considered_spots;
        range_change_ave.at(ibeam) = rshifter_steps.select_range_shifter(range_change_ave.at(ibeam));
        accu_spots += spots_per_field.at(ibeam);
    }

    std::vector<float> remaining_range_shift(nbeams, 0);
    for (size_t i = 0; i < nbeams; i++) {
        if (rs.at(i).exists) {
            // If there is a range add all you need by changing its density.
            // If you need to substract, then stop when you have nothing left and see if you have
            // range remaining.
            if (range_change_ave.at(i) > 0) {
                if (range_change_ave.at(i) <= rs.at(i).wepl) {
                    rs.at(i).substract(range_change_ave.at(i));
                } else {
                    rs.at(i).substract(rs.at(i).wepl);
                    remaining_range_shift.at(i) = range_change_ave.at(i) - rs.at(i).wepl;
                }
            } else {
                rs.at(i).add(range_change_ave.at(i));
            }
        } else {
            // If there is none but one can be provided, create it. Else, add range shift
            // to energies
            if (range_change_ave.at(i) < 0) {
                rs.at(i).create(isocenter_to_beam_distance.at(i), range_change_ave.at(i));
            } else {
                remaining_range_shift.at(i) = range_change_ave.at(i);
            }
        }
    }

    accu_spots = 0;
    for (size_t i = 0; i < nbeams; i++) {        
        if (remaining_range_shift.at(i) != 0) {
            // what is the energy shift per layer so that we get the remaining range shift?
            for (short ispot = 0; ispot < spots_per_field.at(i); ispot++) {
                const float& E0 = orig_energies.at(ispot + accu_spots);
                const float R0 = e_r[0] + e_r[1]*E0 + e_r[2]*std::pow(E0, 2) +
                                 e_r[3]*std::pow(E0, 3) + e_r[4]*std::pow(E0, 4) +
                                 e_r[5]*std::pow(E0, 5);
                const float R = R0 + remaining_range_shift.at(i);
                const float E1 = std::pow(R/r_e[0], 1/r_e[1]) + std::pow(R/r_e[2], 1/r_e[3]) +
                                 std::pow(R/r_e[4], 1/r_e[5]);
                new_energies.at(ispot + accu_spots) = E1*1e6;
            }
            accu_spots += spots_per_field.at(i);
        } else {
            for (short ispot = 0; ispot < spots_per_field.at(i); ispot++) {
                new_energies.at(ispot + accu_spots) = orig_energies.at(ispot + accu_spots);
            }
        }
    }
}


void correct_energy_range (const RShifter_steps_t& rshifter_steps,
                           std::vector<float>& energies,
                           std::vector<RangeShifter_Dims_t>& rs,
                           const std::vector<float>& isocenter_to_beam_distance,
                           const std::vector<short>& spots_per_field,
                           const std::vector<float>& below_thres,
                           const std::string& machine)
{
    const double e_r[6] = {-0.04796, 0.008392, 0.0008097, -1.357e-06, 2.099e-09, -1.754e-12};
    const double r_e[6] = {7.99143832, 6.69077816e-01, 1.96456917e-03,
                           1.80100513, 7.32441249e-08, -2.11977495e-01};
    // Get minimum energy and range for machine
    float min_e{}, max_e{};
    if (machine == "topassmallspots"     ||
        machine == "topas_a5_smallspots" ||
        machine == "topasmediumSpots"    ||
        machine == "topas_a5_mediumspots") {
        min_e = 59.769*1e6;
        max_e = 230.527*1e6;
    } else if (machine == "topasmghr4" ||
               machine == "topasmghr5") {
        min_e = 91.015*1e6;
        max_e = 223.58*1e6;
    }
    float min_range = e_r[0] + e_r[1]*min_e/1e6 + e_r[2]*std::pow(min_e/1e6, 2) +
                      e_r[3]*std::pow(min_e/1e6, 3) + e_r[4]*std::pow(min_e/1e6, 4) +
                      e_r[5]*std::pow(min_e/1e6, 5);
    
    // Get the range difference so that the energy is deliverable
    uint accu_spots{};
    for (size_t ibeam = 0; ibeam < spots_per_field.size(); ibeam++) {
        if (below_thres.at(ibeam) > 0.f) {
            std::vector<float> temp(energies.begin() + accu_spots,
                                    energies.begin() + accu_spots + spots_per_field.at(ibeam));
            const float& E = below_thres.at(ibeam);
            const float R = e_r[0] + e_r[1]*E/1e6 + e_r[2]*std::pow(E/1e6, 2) +
                            e_r[3]*std::pow(E/1e6, 3) + e_r[4]*std::pow(E/1e6, 4) +
                            e_r[5]*std::pow(E/1e6, 5);
            // Discretize range shifter according to user option
            int rshifter_index = -1;
            bool no_errors = false;
            while (rshifter_index != 0 && no_errors == false) {
                no_errors = true;
                float rshifter;
                if (rs.at(ibeam).exists) {
                    rshifter = rshifter_steps.select_range_shifter(min_range-R+rs.at(ibeam).wepl,
                                                                   rshifter_index);
                } else {
                    rshifter = rshifter_steps.select_range_shifter(min_range-R, rshifter_index);
                }
                rshifter_index = rshifter_steps.get_wepl_index(rshifter);
                // Increase spot energies to compensate for the new range shifter
                for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++) {
                    const float& E0 = temp.at(ispot + accu_spots);
                    const float R0 = e_r[0] + e_r[1]*E0/1e6 + e_r[2]*std::pow(E0/1e6, 2) +
                                     e_r[3]*std::pow(E0/1e6, 3) + e_r[4]*std::pow(E0/1e6, 4) +
                                     e_r[5]*std::pow(E0/1e6, 5);
                    const float R = R0 + rshifter;
                    const float E1 = std::pow(R/r_e[0], 1/r_e[1]) + std::pow(R/r_e[2], 1/r_e[3]) +
                                     std::pow(R/r_e[4], 1/r_e[5]);
                    // std::cerr << "Energy of spot @ beam: " << ispot << " @ " << ibeam << " = ";
                    // std::cerr << temp.at(ispot + accu_spots)/1e6 << " MeV ";
                    // std::cerr << "changed to " << E1 << " MeV." << std::endl;
                    if (E1*1e6 > max_e) {
                        std::cerr << "WARNING! The range shifter to force energies too small to ";
                        std::cerr << "be in the gantry energy range makes the high energy spots ";
                        std::cerr << "outside the possible energies: " << E1 << " > " << max_e/1e6;
                        std::cerr << ". The distal edge has higher importance so this is not ";
                        std::cerr << "allowed. A thinner ranger shifter will be attempted, ";
                        std::cerr << "reducing the capping effect of the underflow spots of beam ";
                        std::cerr << "added or increased for beam " << ibeam << "." << std::endl;
                        no_errors = false;
                        break;
                    }
                    temp.at(ispot + accu_spots) = E1*1e6;
                }
                if (no_errors) {
                    std::copy(temp.begin(), temp.end(), energies.begin() + accu_spots);
                    // Create range shifter to make it deliverable
                    std::cerr << "New range shifter for beam " << ibeam << " = ";
                    std::cerr << rshifter << " water cm." << std::endl;
                    if (rs.at(ibeam).exists) {
                        rs.at(ibeam).set_wepl(rshifter);
                    } else {
                        rs.at(ibeam).create(isocenter_to_beam_distance.at(ibeam), rshifter);
                    }
                } else {
                    // There were errors, cap energy to minimum!
                    for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++) {
                        energies.at(ispot+accu_spots) = std::max(energies.at(ispot+accu_spots),
                                                                 min_e);
                    }
                }
            }
        }
        accu_spots += spots_per_field.at(ibeam);
    }
}

