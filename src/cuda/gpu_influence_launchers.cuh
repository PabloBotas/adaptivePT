#ifndef __GPU_INFLUENCE_LAUNERS_CUH__
#define __GPU_INFLUENCE_LAUNERS_CUH__

#include "patient_parameters.hpp"
#include "special_types.hpp"
#include "vector4.hpp"
#include "volume.hpp"

#include <string>
#include <vector>

void influence_from_beam_model_launcher(Array4<float>& influence,
    const std::vector<float>& new_energies,
    const Volume_metadata_t& ct_metadata,
    const Patient_Parameters_t& patient_parameters,
    const uint nspots, const uint n_probing_positions,
    const std::string outputdir);

#endif