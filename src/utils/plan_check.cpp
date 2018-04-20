#include "plan_check.hpp"

#include "volume.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <valarray>
using Array = std::valarray<float>;

// Dose stats -----------------------------------------------
DoseStats::DoseStats (uint n_) : n(n_) {
}
void DoseStats::add_voxel_dose (float dose) {
    min = std::min(min, dose);
    max = std::max(max, dose);
    mean += dose/n;
}
// Dose stats for OARs --------------------------------------
DoseStatsOARs::DoseStatsOARs (uint n_) : DoseStats(n_) {
}
void DoseStatsOARs::check_plan () {
    std::cout << "The adaptation process has the following dose levels in OARs:" << std::endl;
    std::cout << "\tMin dose:  " << min << " Gy" << std::endl;
    std::cout << "\tMean dose: " << mean << " Gy" << std::endl;
    std::cout << "\tMax dose:  " << max << " Gy" << std::endl;
}
// Dose stats for target ------------------------------------
DoseStatsTarget::DoseStatsTarget (uint n_, float prescription_) :
    DoseStats(n_), over(n_), over98(n_), over95(n_), over90(n_),
    over105(n_), over107(n_), over110(n_), over120(n_), prescription(prescription_) {};
void DoseStatsTarget::add_voxel_dose (float dose) {
    min = std::min(min, dose);
    max = std::max(max, dose);
    mean += dose/n;
    if (dose < 0.90*prescription) {
        over90--; over95--; over98--; over--;
        over105--; over107--; over110--; over120--;
    } else if (dose < 0.95*prescription) {
        over95--; over98--; over--;
        over105--; over107--; over110--; over120--;
    } else if (dose < 0.98*prescription) {
        over98--; over--;
        over105--; over107--; over110--; over120--;
    } else if (dose < prescription) {
        over--;
        over105--; over107--; over110--; over120--;
    } else if (dose < 1.05*prescription) {
        over105--; over107--; over110--; over120--;
    } else if (dose < 1.07*prescription) {
        over107--; over110--; over120--;
    } else if (dose < 1.10*prescription) {
        over110--; over120--;
    } else if (dose < 1.20*prescription) {
        over120--;
    }
}
float DoseStatsTarget::check_plan () {
    ratio = float(over)/float(n);
    std::cout << "The adaptation process kept the following coverage:" << std::endl;
    std::cout << "\tMin dose in target:  " << min << " Gy" << std::endl;
    std::cout << "\tMean dose in target: " << mean << " Gy" << std::endl;
    std::cout << "\tMax dose in target:  " << max << " Gy" << std::endl;
    std::cout << "\t" << float(over90)/n*100 << " % above 90% prescription" << std::endl;
    std::cout << "\t" << float(over95)/n*100 << " % above 95% prescription" << std::endl;
    std::cout << "\t" << float(over98)/n*100 << " % above 98% prescription" << std::endl;
    std::cout << "\t" << ratio*100 << " % above prescription" << std::endl;
    std::cout << "\t" << float(over105)/n*100 << " % above 105% prescription" << std::endl;
    std::cout << "\t" << float(over107)/n*100 << " % above 107% prescription" << std::endl;
    std::cout << "\t" << float(over110)/n*100 << " % above 110% prescription" << std::endl;
    std::cout << "\t" << float(over120)/n*100 << " % above 120% prescription" << std::endl;
    return ratio;
}

void check_adaptation_from_dose (std::string adapt_total_dose_file,
                                 const Volume_t& target_mask, const Volume_t& oars_mask,
                                 float dose_prescription, float conv_factor)
{
    std::cout << "Assessing the adaptation!!" << std::endl;
    // Select underdose areas of total dose in target and get underdosage percentage
    auto glambda = [](const float& i){return i > 0.5;};
    uint target_nvox = std::count_if(target_mask.data.begin(), target_mask.data.end(), glambda);
    uint oars_nvox = std::count_if(oars_mask.data.begin(), oars_mask.data.end(), glambda);
    uint volume_vox = target_mask.nElements;

    DoseStatsTarget target_stats(target_nvox, dose_prescription);
    DoseStatsOARs oar_stats(oars_nvox);
    // Get total dose in target
    Volume_t total_adapt_dose(adapt_total_dose_file, Volume_t::Source_type::DOSE);
    total_adapt_dose.scale(conv_factor);
    for (size_t i = 0; i < volume_vox; ++i) {
        if (target_mask.data.at(i) > 0.5) {
            target_stats.add_voxel_dose(total_adapt_dose.data[i]);
        }
        if (oars_mask.data.at(i) > 0.5) {
            oar_stats.add_voxel_dose(total_adapt_dose.data[i]);
        }
    }

    float res = target_stats.check_plan();
    oar_stats.check_plan();
    if (res > 0.98) {
        std::cout << "This is great!! Shooooooot!!" << std::endl;
    }
}


void check_adaptation (const Array& dose,
                       float dose_prescription,
                       uint target_nvox, uint oars_nvox, uint dij_nvox,
                       const std::vector<float>& total_mask,
                       const Volume_t& target_mask,
                       const Volume_t& oars_mask,
                       std::string out_directory)
{
    // Get stats!!
    DoseStatsTarget target_stats(target_nvox, dose_prescription);
    DoseStatsOARs oar_stats(oars_nvox);
    for (size_t i = 0, imask = 0; i < total_mask.size() && imask < dij_nvox; ++i) {
        // Check if in a ROI
        if (total_mask.at(i) > 0.5) {
            if (target_mask.data.at(i) > 0.5) {
                target_stats.add_voxel_dose(dose[i]);
            }
            if (oars_mask.data.at(i) > 0.5) {
                oar_stats.add_voxel_dose(dose[i]);
            }
            imask++;
        }
    }

    std::ofstream ofs(out_directory+"/opt4D_optimized_dose.dat", std::ios::binary);
    ofs.write((char*)&dose[0], dose.size()*sizeof(float));
    ofs.close();

    std::cout << std::endl;
    float res = target_stats.check_plan();
    oar_stats.check_plan();
    if (res > 0.98) {
        std::cout << "This is great!! Shooooooot!!" << std::endl;
    }
}
