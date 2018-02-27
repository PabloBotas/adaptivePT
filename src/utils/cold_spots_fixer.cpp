#include "cold_spots_fixer.hpp"

#include "opt4D_manager.hpp"
#include "utils.hpp"
#include "volume.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <utility>
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
                       Array& dose_in_target)
{
    std::cout << "Assessing the adaptation!!" << std::endl;
    // Select underdose areas of total dose in target and get underdosage percentage
    auto glambda = [](const float& i){return i > 0.5;};
    uint target_nvox = std::count_if(target_mask.data.begin(), target_mask.data.end(), glambda);
    // uint target_rim_nvox = std::count_if(target_rim_mask.data.begin(),
                                        // target_rim_mask.data.end(), glambda);
    uint oars_nvox = std::count_if(oars_mask.data.begin(), oars_mask.data.end(), glambda);
    // There are overlaps, so the total size is not the sum of the individuals
    std::vector<float> total_mask = utils::or_vectors(target_mask.data, target_rim_mask.data,
                                                      oars_mask.data);
    uint dij_nvox = std::count_if(total_mask.begin(), total_mask.end(), glambda);
    uint volume_vox = target_mask.nElements;
    uint nbeams = adapt_dij_per_field_files.size();

    // Fill total dose, dose per beam, underdose_mask, dose in target
    dose.resize(volume_vox, 0);
    adapt_field_dose.resize(nbeams);
    for (size_t i = 0; i < nbeams; ++i)
        adapt_field_dose.at(i).resize(volume_vox, 0);
    dose_in_target.resize(target_nvox);
    underdose_mask.resize(volume_vox);
    dose_in_mask.resize(dij_nvox);

    // Read dijs inside masked region
    field_dij.resize(nbeams);
    std::cout << "Read dijs inside target:" << std::endl;
    for (size_t i = 0; i < nbeams; ++i) {
        // Read Dij file
        std::string& f = adapt_dij_per_field_files.at(i);
        std::cout << "\t- " << f << std::endl;
        std::ifstream ifs(f, std::ios::binary);
        utils::check_fs(ifs, f, " to read Dij during plan assessment");
        field_dij.at(i).resize(spots_per_field.at(i)*dij_nvox);
        ifs.read(reinterpret_cast<char*>(&field_dij.at(i)[0]),
                                         spots_per_field.at(i)*dij_nvox*sizeof(float));
        field_dij.at(i) *= conv_factor;
    }

    DoseStatsTarget target_stats(target_nvox, dose_prescription);
    DoseStatsOARs oar_stats(oars_nvox);
    // Iterate over volume voxels
    for (size_t i = 0, idij = 0, itarget = 0; i < volume_vox; ++i) {
        // Check if in a ROI
        if (total_mask.at(i) > 0.5) {
            for (size_t f = 0; f < nbeams; ++f) {
                // Get dose in voxel per spot
                Array d_per_spot = field_dij.at(f)[
                        std::slice(idij, spots_per_field.at(f), dij_nvox)];
                float d = d_per_spot.sum();
                dose[i] += d;
                adapt_field_dose.at(f)[i] += d;
                dose_in_mask[idij] += d;
            }
            idij++;
            if (target_mask.data.at(i) > 0.5) {
                dose_in_target[itarget] = dose[i];
                if (dose[i] < dose_prescription)
                    underdose_mask[i] = dose_prescription - dose[i];
                target_stats.add_voxel_dose(dose[i]);
                itarget++;
            }
            if (oars_mask.data.at(i) > 0.5) {
                oar_stats.add_voxel_dose(dose[i]);
            }
        }
    }

    // for (size_t ibeam = 0; ibeam < nbeams; ++ibeam) {
    //     for (short ispot = 0; ispot < spots_per_field.at(ibeam); ++ispot) {
    //         Array d = field_dij.at(ibeam)[
    //                                                 std::slice(ispot*target_nvox, target_nvox, 1)];
    //         std::cout << "Average dose beam " << ibeam << " spot " << ispot;
    //         std::cout << std::setprecision(10) << " = " << d.sum()/target_nvox << std::endl;
    //     }
    // }

    float res = target_stats.check_plan();
    oar_stats.check_plan();
    if (res > 0.98) {
        std::cout << "This is great!! Shooooooot!!" << std::endl;
    }
}


void check_adaptation(const std::vector<Array>& adapt_field_dij,
                      float dose_prescription, std::vector<short> spots_per_field,
                      uint target_nvox, uint oars_nvox, uint dij_nvox,
                      const std::vector<float>& total_mask,
                      const Volume_t& target_mask,
                      const Volume_t& oars_mask,
                      std::string out_directory)
{
    // Get total dose from new Dij
    Array adapt_dose_mask(dij_nvox);
    for (size_t f = 0; f < adapt_field_dij.size(); ++f) {
        for (size_t i = 0; i < dij_nvox; ++i) {
            Array d_to_vox_per_spot = adapt_field_dij.at(f)[
                                    std::slice(i, spots_per_field.at(f), dij_nvox)];
            adapt_dose_mask[i] += d_to_vox_per_spot.sum();
        }
    }

    // Get stats!!
    DoseStatsTarget target_stats(target_nvox, dose_prescription);
    DoseStatsOARs oar_stats(oars_nvox);
    Array out;
    out.resize(target_mask.nElements);
    for (size_t i = 0, imask = 0; i < total_mask.size(); ++i) {
        // Check if in a ROI
        if (total_mask.at(i) > 0.5) {
            out[i] = adapt_dose_mask[imask];
            if (target_mask.data.at(i) > 0.5) {
                target_stats.add_voxel_dose(adapt_dose_mask[imask]);
            }
            if (oars_mask.data.at(i) > 0.5) {
                oar_stats.add_voxel_dose(adapt_dose_mask[imask]);
            }
            imask++;
        }
    }

    std::ofstream ofs(out_directory+"/opt4D_optimized_dose.dat", std::ios::binary);
    ofs.write((char*)&out[0], out.size()*sizeof(float));
    ofs.close();

    std::cout << std::endl;
    float res = target_stats.check_plan();
    oar_stats.check_plan();
    if (res > 0.98) {
        std::cout << "This is great!! Shooooooot!!" << std::endl;
    }
}


void cold_spots_fixer(const Array& adapt_dose_in_mask,           // total dose
                      std::vector<Array>& adapt_field_dij,       // dose per spot
                      const Volume_t& target_mask,               // target region
                      const Volume_t& rim_mask,                  // target rim region
                      const Volume_t& oars_mask,                 // OARs region
                      float dose_prescription,                   // prescription
                      std::string out_directory,                 // output dir
                      const std::vector<short>& spots_per_field, // # of spots per field
                      std::vector<float> weights,                // spots weights
                      std::vector<float>& weight_scaling)        // spots scaling (all 1)
{
    // uint volume_vox = target_mask.nElements;
    uint nbeams = adapt_field_dij.size();
    uint total_spots = weights.size();
    std::vector<uint> accu_spots(spots_per_field.size());
    std::partial_sum(spots_per_field.begin(), spots_per_field.end(), accu_spots.begin());

    auto glambda = [](const float& i) {return i > 0.5;};
    uint target_nvox = std::count_if(target_mask.data.begin(), target_mask.data.end(), glambda);
    uint oars_nvox = std::count_if(oars_mask.data.begin(), oars_mask.data.end(), glambda);
    uint rim_nvox = std::count_if(rim_mask.data.begin(), rim_mask.data.end(), glambda);

    std::vector<float> total_mask = utils::or_vectors(target_mask.data, rim_mask.data,
                                                      oars_mask.data);
    uint dij_nvox = std::count_if(total_mask.begin(), total_mask.end(), glambda);

    // Select spots for weight adjustment!!
    // The spots must be well represented in the Dij. I could read the std file and deduce
    // from there, but that would force me to go vox by vox and do a coparison. Let's try
    // first with the weight. Let's sort all the spots by weight and select the higher
    // weighted ones with some statistical method. Which?
    //  - select the spots accounting for 25% of the total weight? 50%? 75%?
    //  - From patient 3: 20% of the spots account for 65% of the dose and 58% of the weight
    //  - Cummulative dose and weight are essentially the same curve for our purposes
    //  - 50% of the weight corresponds to 15% of the spots. I will start with that, rounding up

    // 1: Create array with weight, index and sort it. Sort will take the first item of the pair
    std::vector<std::pair<float, uint>> indexed_w;
    indexed_w.reserve(total_spots);
    for (uint i = 0; i < weights.size(); ++i)
        indexed_w.push_back(std::make_pair(weights.at(i), i));
    std::sort(indexed_w.begin(), indexed_w.end(), std::greater<std::pair<float, uint>>());
    std::sort(weights.begin(), weights.end(), std::greater<float>());

    float total_w = std::accumulate(weights.begin(), weights.end(), 0.0);
    // 2: Select spots giving 50% of the weight (Cummulative sum IN PLACE!!)
    std::partial_sum(weights.begin(), weights.end(), weights.begin());
    const float reopt_weight_fraction = 0.5;
    uint selected_spots = std::distance(weights.begin(),
                                        std::upper_bound(weights.begin(), weights.end(),
                                        reopt_weight_fraction*total_w));
    std::vector<uint> indexes(selected_spots);
    for (uint i = 0; i < selected_spots; ++i)
        indexes.at(i) = indexed_w.at(i).second;
    // It's convenient to sort them again for clarity
    std::sort(indexes.begin(), indexes.end());
    // Count number of selected spots per field
    std::vector<uint> selected_spots_per_field(nbeams);
    for (size_t i = 0; i < indexes.size(); ++i) {
        uint ibeam = std::distance(accu_spots.begin(),
                                   std::upper_bound(accu_spots.begin(), accu_spots.end(),
                                                    indexes.at(i)));
        selected_spots_per_field.at(ibeam)++;
    }

    std::cout << "Considering " << selected_spots << " spots (";
    std::cout << 100*selected_spots/float(total_spots) << " %), ";
    std::cout << 100*weights.at(selected_spots)/total_w << " % of the weight. ";
    std::cout << "There is an almost linear relation between spot weight and dose." << std::endl;

    // 3: Subset Dij to pass to Opt4D
    // Allocate size in advance. Resizing a valarray deallocates all the items in it, erasing them
    std::cout << "Preparing to subset Dij..." << std::endl;
    std::vector<Array> subset_dij(nbeams);
    for (size_t i = 0; i < nbeams; ++i) {
        subset_dij.at(i).resize(selected_spots_per_field.at(i)*dij_nvox);
    }
    // Subset
    std::cout << "Subsetting Dij..." << std::endl;
    for (size_t i = 0, j_spot_beam = 0; i < indexes.size(); ++i) {
        uint ibeam = std::distance(accu_spots.begin(),
                                   std::upper_bound(accu_spots.begin(), accu_spots.end(),
                                                    indexes.at(i)));
        uint ispot = ibeam > 0 ? indexes.at(i) - accu_spots.at(ibeam-1) : indexes.at(i);
        // Get start position of slice and add it!!
        subset_dij.at(ibeam)[std::slice(j_spot_beam*dij_nvox, dij_nvox, 1)] = 
                                adapt_field_dij.at(ibeam)[std::slice(ispot*dij_nvox, dij_nvox, 1)];
        j_spot_beam = j_spot_beam == selected_spots_per_field.at(ibeam)-1 ? 0 : j_spot_beam+1;
    }

    // 4: Create target dose distribution
    // The target dose is the prescribed dose minus the dose by the spots not-included for
    // reoptimization. So, I substract the dose by the included spots to the total and then take
    // that from the prescription. That will be the dose the spots subset have to provide.
    std::cout << "Calculating target dose..." << std::endl;
    Array target_dose = adapt_dose_in_mask;
    float zero = 0;
    for (size_t i = 0, idij = 0;
         i < target_mask.nElements && idij < dij_nvox;
         ++i) {
        if (total_mask.at(i) > 0.5) {
            if (target_mask.data.at(i) > 0.5) {
                for (size_t f = 0; f < nbeams; ++f) {
                    Array temp = subset_dij[f][
                        std::slice(idij, selected_spots_per_field.at(f), dij_nvox)];
                    target_dose[idij] -= temp.sum();
                    target_dose[idij] = std::max(zero, target_dose[idij]);
                }
                target_dose[idij] = dose_prescription - target_dose[idij];
            }
            idij++;
        }
    }

    // Array target_dose_vol;
    // target_dose_vol.resize(target_mask.nElements);
    // for (size_t i = 0, idij = 0;
    //      i < target_mask.nElements && idij < dij_nvox;
    //      ++i) {
    //     if (total_mask.at(i) > 0.5) {
    //         if (target_mask.data.at(i) > 0.5) {
    //             target_dose_vol[i] = target_dose[idij];
    //         }
    //         idij++;
    //     }
    // }
    // std::ofstream ofs(out_directory+"/dose_mask_opt4D.dat", std::ios::binary);
    // ofs.write((char*)&target_dose[0], target_dose.size()*sizeof(float));
    // ofs.close();

    // 5: Prepare and launch Opt4D and get the weigth scaling
    Opt4D_manager opt4d(out_directory);
    opt4d.populate_directory(selected_spots, dij_nvox, target_nvox, rim_nvox, oars_nvox,
                             total_mask, target_mask, rim_mask, oars_mask,
                             target_dose, subset_dij);
    opt4d.launch_optimization();
    std::vector<float> subset_weight_scaling = opt4d.get_weight_scaling();

    // 6: Apply weight scaling and set in referenced vector
    for (size_t i = 0; i < indexes.size(); ++i) {
        uint ibeam = std::distance(accu_spots.begin(),
                                   std::upper_bound(accu_spots.begin(), accu_spots.end(),
                                                    indexes.at(i)));
        uint ispot = ibeam > 0 ? indexes.at(i) - accu_spots.at(ibeam-1) : indexes.at(i);
        // Get start position of slice and add it!!
        Array temp = adapt_field_dij.at(ibeam)[std::slice(ispot*dij_nvox, dij_nvox, 1)];
        temp *= subset_weight_scaling.at(i);
        adapt_field_dij.at(ibeam)[std::slice(ispot*dij_nvox, dij_nvox, 1)] = temp;
        // Set in referenced vector
        weight_scaling.at(indexes.at(i)) = subset_weight_scaling.at(i);
    }

    // 7: Check new plan dose
    check_adaptation(adapt_field_dij, dose_prescription, spots_per_field,
                     target_nvox, oars_nvox, dij_nvox, total_mask, target_mask, oars_mask,
                     out_directory);
}


void output_debug_doses (std::string out_directory, const Array& underdose_mask,
                         const Volume_t& target_mask,
                         const std::vector<std::string>& plan_dose_per_field_files,
                         const Array& adapt_dose,
                         const std::vector<Array>& adapt_field_dose,
                         const Array& adapt_dose_in_target,
                         float dose_prescription)
{
    std::ofstream ofs0(out_directory+"/adapt_dose.dat", std::ios::binary);
    ofs0.write((char*)&adapt_dose[0], adapt_dose.size()*sizeof(float));
    ofs0.close();

    Array temp;
    temp.resize(adapt_dose.size());
    for (size_t i = 0, itarget = 0; i < adapt_dose.size(); ++i) {
        if (target_mask.data.at(i) > 0.5) {
            temp[i] = adapt_dose_in_target[itarget];
            itarget++;
        }
    }
    // std::ofstream ofs5(out_directory+"/adapt_dose_to_target.dat", std::ios::binary);
    // ofs5.write((char*)&temp[0], temp.size()*sizeof(float));
    // ofs5.close();

    std::ofstream ofs(out_directory+"/cold_spots_map.dat", std::ios::binary);
    ofs.write((char*)&underdose_mask[0], underdose_mask.size()*sizeof(float));
    ofs.close();

    // std::ofstream ofs2(out_directory+"/target_mask.dat", std::ios::binary);
    // ofs2.write((char*)target_mask.data.data(), target_mask.data.size()*sizeof(float));
    // ofs2.close();

    // uint nbeams = plan_dose_per_field_files.size();
    // uint volume_vox = target_mask.nElements;
    // std::vector<std::vector<float>> diffs(nbeams, std::vector<float>(volume_vox, 0));
    // Array total_plan_dose;
    // total_plan_dose.resize(volume_vox);
    // for (size_t i = 0; i < nbeams; ++i) {
    //     Volume_t plan_dose(plan_dose_per_field_files.at(i), Volume_t::Source_type::DOSE);
    //     for (size_t j = 0; j < volume_vox; ++j) {
    //         total_plan_dose[j] += plan_dose.data.at(j);
    //         if (target_mask.data.at(j) > 0.5 &&
    //             adapt_dose[j] < dose_prescription &&
    //             adapt_field_dose.at(i)[j] < plan_dose.data.at(j)) {
    //             // Fill diff per beam with cold spot info
    //             diffs.at(i).at(j) = plan_dose.data.at(j) - adapt_field_dose.at(i)[j];
    //         }
    //     }
    //     std::ofstream ofs3(out_directory+"/cold_spots_map.beam_" +
    //                        std::to_string(i)+".dat", std::ios::binary);
    //     ofs3.write((char*)diffs.at(i).data(), diffs.at(i).size()*sizeof(float));
    //     ofs3.close();
    // }
}
