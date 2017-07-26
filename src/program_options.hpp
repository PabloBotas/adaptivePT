#ifndef __WARP_OPTIONS__
#define __WARP_OPTIONS__

#include "special_types.hpp"
#include <vector>

enum Warp_opts_t
{
    FREE_POS_FREE_ENERGY,
    FREE_POS_RIGID_ENERGY,
    FREE_POS_RIGID_BEAMS_ENERGY,

    RIGID_POS_FREE_ENERGY,
    RIGID_POS_RIGID_ENERGY,
    RIGID_POS_RIGID_BEAMS_ENERGY,

    RIGID_BEAMS_POS_FREE_ENERGY,
    RIGID_BEAMS_POS_RIGID_ENERGY,
    RIGID_BEAMS_POS_RIGID_BEAMS_ENERGY,

    // RIGID_POS_LAYERS_FREE_ENERGY,
    // RIGID_ENERGY_LAYERS_FREE_POS,
    // RIGID_ENERGY_LAYERS_RIGID_POS_LAYERS,
};

void apply_energy_options (Warp_opts_t opts,
                           std::vector<double>& energy_shift,
                           const std::vector<short>& spots_per_field);
void apply_rigid_energy (std::vector<double>& energy_shift);
void apply_rigid_energy_per_beam (std::vector<double>& energy,
                                  const std::vector<short>& spots_per_field);
void apply_position_options (Warp_opts_t opts,
                             Array3<double>& vf,
                             const std::vector<short>& spots_per_field);
void apply_rigid_positions (Array3<double>& vf);
void apply_rigid_positions_per_beam (Array3<double>& vf,
                                     const std::vector<short>& spots_per_field);

#endif