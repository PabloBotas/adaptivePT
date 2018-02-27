#ifndef __COLD_SPOTS_FIXER_HPP__
#define __COLD_SPOTS_FIXER_HPP__

#include "utils.hpp"
#include "volume.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <valarray>
using Array = std::valarray<float>;

struct DoseStats {
    uint n;
    float min = 1E7;
    float max{};
    float mean{};
    DoseStats (uint n_);
    void add_voxel_dose (float dose);
};

struct DoseStatsOARs : DoseStats {
    DoseStatsOARs (uint n_);
    void check_plan ();
};

struct DoseStatsTarget : DoseStats {
    uint over;
    uint over98;
    uint over95;
    uint over90;
    uint over105;
    uint over107;
    uint over110;
    uint over120;
    float ratio;
    float prescription;
    DoseStatsTarget (uint n_, float prescription_);
    void add_voxel_dose (float dose);
    float check_plan ();
};

void check_adaptation_from_dose (std::string adapt_total_dose_file,
                                 const Volume_t& target_mask, const Volume_t& oars_mask,
                                 float dose_prescription, float conv_factor);

void check_adaptation (std::vector<std::string> adapt_dij_per_field_files,
                       float dose_prescription, float conv_factor,
                       const std::vector<short>& spots_per_field,
                       std::vector<Array>& field_dij,
                       Array& dose,
                       std::vector<Array>& adapt_field_dose,
                       Array& underdose_mask,                      
                       const Volume_t& target_mask,
                       const Volume_t& target_rim_mask,
                       const Volume_t& oars_mask,
                       Array& dose_in_mask,
                       Array& dose_in_target);

void check_adaptation(const std::vector<Array>& adapt_field_dij,
                      float dose_prescription, std::vector<short> spots_per_field,
                      uint target_nvox, uint oars_nvox, uint dij_nvox,
                      const std::vector<float>& total_mask,
                      const Volume_t& target_mask,
                      const Volume_t& oars_mask,
                      std::string out_directory);

void cold_spots_fixer(const Array& adapt_dose_in_mask,           // total dose
                      std::vector<Array>& adapt_field_dij,       // dose per spot
                      const Volume_t& target_mask,               // target region
                      const Volume_t& rim_mask,                  // target rim region
                      const Volume_t& oars_mask,                 // OARs region
                      float dose_prescription,                   // prescription
                      std::string out_directory,                 // output dir
                      const std::vector<short>& spots_per_field, // # of spots per field
                      std::vector<float> weights,                // spots weights
                      std::vector<float>& weight_scaling);       // spots scaling (all 1)

void output_debug_doses (std::string out_directory, const Array& underdose_mask,
                         const Volume_t& target_mask,
                         const std::vector<std::string>& plan_dose_per_field_files,
                         const Array& adapt_dose,
                         const std::vector<Array>& adapt_field_dose,
                         const Array& adapt_dose_in_target,
                         float dose_prescription);

#endif