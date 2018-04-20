#ifndef __PLAN_CHECK_HPP__
#define __PLAN_CHECK_HPP__

#include "volume.hpp"

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

void check_adaptation (const Array& dose,
                       float dose_prescription,
                       uint target_nvox, uint oars_nvox, uint dij_nvox,
                       const std::vector<float>& total_mask,
                       const Volume_t& target_mask,
                       const Volume_t& oars_mask,
                       std::string out_directory);

#endif