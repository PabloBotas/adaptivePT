#ifndef __WEIGHT_ADAPTATION_HPP__
#define __WEIGHT_ADAPTATION_HPP__

#include "volume.hpp"

#include <string>
#include <vector>
#include <valarray>
using Array = std::valarray<float>;

void select_spots_by_weight (uint& nspots_subset,
                             std::vector<std::vector<uint>>& indices,
                             std::vector<uint>& subset_per_field,
                             const std::vector<float>& weights,
                             const std::vector<uint>& accu_spots,
                             const uint nbeams,
                             const float weight_frac = 0.5,
                             const float spots_frac = 0.1);

Array read_transpose_dij (const std::string& file,
                          const uint nspots,
                          const uint nvox,
                          const float conv_factor);

void fill_dose_arrays (const Array& dij, const uint dij_nvox,
                       const uint nspots, const uint volume_vox,
                       const std::vector<float>& mask, const Volume_t& target,
                       Array& dose, Array& dose_in_target, Array& dose_in_mask);

void adapt_weights (const std::vector<std::string>& dij_files,         // Dij files
                    const std::vector<std::string>& target_mask_files, // target region
                    const std::vector<std::string>& rim_mask_files,    // target rim region
                    const std::vector<std::string>& oars_mask_files,   // OARs region
                    const float dose_prescription,                     // prescription
                    const float conv_factor,                           // factor to Gy
                    const std::string& out_directory,                  // output dir
                    const std::vector<short>& spots_per_field,         // # of spots per field
                    std::vector<float> weights,                        // spots weights
                    std::vector<float>& weight_scaling);               // spots scaling (all 1)

// void output_debug_doses (std::string out_directory, const Array& underdose_mask,
//                          const Volume_t& target_mask,
//                          const std::vector<std::string>& plan_dose_per_field_files,
//                          const Array& adapt_dose,
//                          const std::vector<Array>& adapt_field_dose,
//                          const Array& adapt_dose_in_target,
//                          float dose_prescription);

#endif