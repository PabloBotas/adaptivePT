#include "weight_adaptation.hpp"

#include "plan_check.hpp"
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
#include <omp.h>
#include <time.h>
using Array = std::valarray<float>;


void select_spots_by_weight (uint& nspots_subset,
                             std::vector<std::vector<uint>>& indices,
                             std::vector<uint>& subset_per_field,
                             const std::vector<float>& weights,
                             const std::vector<uint>& accu_spots,
                             const uint nbeams,
                             const float weight_frac,
                             const float spots_frac)
{
    uint total_spots = weights.size();
    // 1: create paired sorted list
    std::vector<std::pair<float, uint>> indexed_w;
    indexed_w.reserve(total_spots);
    for (uint i = 0; i < total_spots; ++i) {
        indexed_w.push_back(std::make_pair(weights.at(i), i));
    }
    // Sort lists. Sorting a pair sorts taking the first field into account
    std::sort(indexed_w.begin(), indexed_w.end(), std::greater<std::pair<float, uint>>());
    std::vector<float> mod_weights = weights;
    std::sort(mod_weights.begin(), mod_weights.end(), std::greater<float>());

    float total_w = std::accumulate(mod_weights.begin(), mod_weights.end(), 0.0);
    // 2: Select spots giving 50% of the weight
    std::partial_sum(mod_weights.begin(), mod_weights.end(), mod_weights.begin());
    uint selected_w = std::distance(mod_weights.begin(),
                                    std::upper_bound(mod_weights.begin(), mod_weights.end(),
                                    weight_frac*total_w));
    uint selected_min = uint(spots_frac*total_spots + 0.5);
    nspots_subset = std::max(selected_w, selected_min);
    // 3: Fill selected indices
    indices.resize(nbeams);
    subset_per_field.resize(nbeams);
    for (uint i = 0; i < nspots_subset; ++i) {
        uint abs_spot = indexed_w.at(i).second;
        uint ibeam = std::distance(accu_spots.begin(),
                                    std::upper_bound(accu_spots.begin(), accu_spots.end(),
                                    abs_spot));
        uint ispot;
        if (ibeam == 0)
            ispot = abs_spot;
        else
            ispot = abs_spot - accu_spots.at(ibeam-1);
        indices.at(ibeam).push_back(ispot);
        subset_per_field.at(ibeam)++;
    }
    // Sort to have the spots grouped by beam ids
    for (uint i = 0; i < nbeams; ++i) {
        std::sort(indices.at(i).begin(), indices.at(i).end());
    }

    std::cout << "Adjusting: " << std::endl;
    std::cout << "    Selected # - Spot frac - Weight frac" << std::endl;
    std::cout << "    Final " << nspots_subset << "\t" << 100*nspots_subset/float(total_spots);
    std::cout << "\t" << 100*mod_weights.at(nspots_subset)/total_w << std::endl;
    std::cout << "    Min.  " << selected_min << "\t" << 100*selected_min/float(total_spots);
    std::cout << "\t" << 100*mod_weights.at(selected_min)/total_w << std::endl;
    std::cout << "    Prev. " << selected_w << "\t" << 100*selected_w/float(total_spots);
    std::cout << "\t" << 100*mod_weights.at(selected_w)/total_w << std::endl;
    std::cout << "    Spots per field: ";
    for (uint i = 0; i < indices.size(); ++i) {
        std::cout << " " << indices.at(i).size();
    }
    std::cout << std::endl;

    std::cout << "There is an almost linear relation between spot weight and dose." << std::endl;
}

Array read_dij (const std::string& file,
                const uint nspots,
                const uint nvox,
                const float conv_factor)
{
    std::ifstream ifs(file, std::ios::binary);
    utils::check_fs(ifs, file, " to read Dij during plan assessment");
    Array dij(nspots*nvox);
    ifs.read(reinterpret_cast<char*>(&dij[0]), nspots*nvox*sizeof(float));
    dij *= conv_factor;
    
    return dij;
}

Array read_transpose_dij (const std::string& file,
                          const uint nspots,
                          const uint nvox,
                          const float conv_factor)
{
    Array temp = read_dij(file, nspots, nvox, conv_factor);
    Array dij(nspots*nvox);
    #pragma omp parallel for
    for (size_t j = 0; j < nvox; ++j) {
        dij[std::slice(j*nspots, nspots, 1)] = temp[std::slice(j, nspots, nvox)];
    }
    return dij;
}

void fill_dose_arrays (const Array& dij, const uint dij_nvox,
                       const uint nspots, const uint volume_vox,
                       const std::vector<float>& mask, const Volume_t& target,
                       Array& dose, Array& dose_in_target, Array& dose_in_mask)
{
    // std::cout << "Creating map!" << std::endl;
    std::vector<int> dij_mask_map(volume_vox, -1);
    for (size_t i = 0, idij = 0;
        i < volume_vox && idij < dij_nvox;
        ++i) {
        if (mask.at(i) > 0.5) {
            dij_mask_map[i] = idij++;
        }
    }

    // uint tot = 0;
    // std::cout << "Subsetting!" << std::endl;
    for (size_t i = 0; i < volume_vox; ++i) {
        // Check if in a ROI
        if (dij_mask_map[i] > -1) {
            // Accumulate dose in voxel per spot
            Array temp = dij[std::slice(dij_mask_map[i], nspots, dij_nvox)];
            float d = temp.sum();
            dose[i] += d;
            dose_in_mask[dij_mask_map[i]] += d;
            // if (target.data.at(i) > 0.5) {
            //     dose_in_target[itarget] += d;
            //     itarget++;
            // }
            // tot++;
        }
    }
    // std::cout << "Total vox " << tot << std::endl;
}

// void fill_dose_arrays (const Array& dij, const uint dij_nvox,
//                        const uint nspots, const uint volume_vox,
//                        const std::vector<float>& mask, const Volume_t& target,
//                        Array& dose, Array& dose_in_target, Array& dose_in_mask)
// {
//     for (size_t i = 0, idij = 0, itarget = 0; i < volume_vox; ++i) {
//         // Check if in a ROI
//         if (mask.at(i) > 0.5) {
//             // Get dose in voxel per spot
//             // auto p = std::begin(dij)+idij*nspots; 
//             // float d = std::accumulate(p, p+nspots, 0.0);
//             Array temp = dij[std::slice(idij, nspots, dij_nvox)];
//             float d = temp.sum();
//             dose[i] += d;
//             dose_in_mask[idij] += d;
//             idij++;
//             if (target.data.at(i) > 0.5) {
//                 dose_in_target[itarget] += d;
//                 itarget++;
//             }
//         }
//     }
// }

void adapt_weights (const std::vector<std::string>& dij_files,         // Dij files
                    const std::vector<std::string>& target_mask_files, // target region
                    const std::vector<std::string>& rim_mask_files,    // target rim region
                    const std::vector<std::string>& oars_mask_files,   // OARs region
                    const float dose_prescription,                     // prescription
                    const float conv_factor,                           // factor to Gy
                    const std::string& out_directory,                  // output dir
                    const std::vector<short>& spots_per_field,         // # of spots per field
                    std::vector<float> weights,                        // spots weights
                    std::vector<float>& weight_scaling)                // spots scaling (all 1)
{
    // 1: Definitions --------------------------------------------------------
    uint nbeams = dij_files.size();
    // uint total_spots = weights.size();
    std::vector<uint> accu_spots(spots_per_field.size());
    std::partial_sum(spots_per_field.begin(), spots_per_field.end(), accu_spots.begin());

    // 2: Select spots -------------------------------------------------------
    uint nspots_subset;
    std::vector<std::vector<uint>> indexes;
    std::vector<uint> subset_per_field;
    select_spots_by_weight(nspots_subset, indexes, subset_per_field, weights, accu_spots, nbeams);

    // 3: Read masks and get data --------------------------------------------
    Volume_t target_mask = utils::read_masks (target_mask_files);
    Volume_t target_rim_mask = utils::read_masks (rim_mask_files);
    Volume_t oars_mask = utils::read_masks (oars_mask_files);
    uint target_nvox = target_mask.count_above_thres(0.5);
    uint rim_nvox = target_rim_mask.count_above_thres(0.5);
    uint oars_nvox = oars_mask.count_above_thres(0.5);
    // There are overlaps, so the total size is not the sum of the individuals
    std::vector<float> total_mask = utils::or_vectors(target_mask.data, target_rim_mask.data,
                                                      oars_mask.data);
    auto glambda = [](const float& i){return i > 0.5;};
    uint dij_nvox = std::count_if(total_mask.begin(), total_mask.end(), glambda);
    uint volume_vox = target_mask.nElements;


    // 4: Read Dijs: get total dose and spots subset -------------------------
    // Fill total dose, dose per beam, dose in target
    Array dose;
    Array dose_in_target;
    Array dose_in_mask;
    dose.resize(volume_vox, 0);
    dose_in_target.resize(target_nvox, 0);
    dose_in_mask.resize(dij_nvox, 0);

    std::vector<Array> subset_dij(nbeams);
    std::cout << "Reading and subsetting dijs ..." << std::endl;
    for (size_t ibeam = 0; ibeam < nbeams; ++ibeam) {
        // Read Dij file
        std::cout << "\t- " << dij_files.at(ibeam);
        Array dij = read_dij(dij_files.at(ibeam), spots_per_field.at(ibeam),
                             dij_nvox, conv_factor);
        std::cout << " ... read!" << std::endl;
        // Get selected spots
        subset_dij.at(ibeam).resize(indexes.at(ibeam).size()*dij_nvox);
        #pragma omp parallel for
        for (size_t i = 0; i < indexes.at(ibeam).size(); ++i) {
            uint& ispot = indexes.at(ibeam).at(i);
            // Copy spot to subset
            std::copy(std::begin(dij) + dij_nvox*ispot,
                      std::begin(dij) + dij_nvox*(ispot+1),
                      std::begin(subset_dij.at(ibeam)) + i*dij_nvox);
        }

        // Iterate over volume voxels and fill dose arrays
        fill_dose_arrays(dij, dij_nvox, spots_per_field.at(ibeam), volume_vox,
                         total_mask, target_mask,
                         dose, dose_in_target, dose_in_mask);
    }

    // 5: Test geometrical plan quality --------------------------------------
    check_adaptation(dose, dose_prescription,
                     target_nvox, oars_nvox, dij_nvox,
                     total_mask, target_mask, oars_mask,
                     out_directory);

    // 6: Create target dose distribution ------------------------------------
    // The target dose is the prescribed dose minus the dose by the spots not-included for
    // reoptimization. So, I substract the dose by the included spots to the total and then take
    // that from the prescription. That will be the dose the spot subset has to provide.
    // The dose missing MUST have the same nvoxels as the Dij to keep correspondence within Opt4D
    // between voxels. The dose missing has zeros outside the target voxels
    std::cout << "Calculating dose missing..." << std::endl;
    Array baseline_dose = dose;
    Array target_dose_missing;
    target_dose_missing.resize(dij_nvox);
    Array dose_by_subset;
    dose_by_subset.resize(volume_vox);
    for (size_t i = 0, idij = 0; i < volume_vox && idij < dij_nvox; ++i) {
        if (total_mask.at(i) > 0.5) {
            float d = 0;
            for (size_t f = 0; f < nbeams; ++f) {
                Array temp = subset_dij[f][std::slice(idij, indexes.at(f).size(), dij_nvox)];
                d += temp.sum();
            }
            dose_by_subset[i] += d;
            baseline_dose[i] -= d;
            if (target_mask.data.at(i) > 0.5) {
                // the result could be negative!! But that would be fine. Negative voxels are set to
                // zero by Opt4D. The reason not to do it here is that I have allowed Opt4D to add a
                // constant value to the voxel-dose to control the maximum dose within the target.
                // This operation needs to know if the voxel is -3 or -10 or whatever value. If
                // after summing, the value is still negative, it will be set to zero.
                target_dose_missing[idij] = dose_prescription - baseline_dose[i];
            }
            idij++;
        }
    }

    // Output debug volumes
    Array target_dose_vol;
    target_dose_vol.resize(volume_vox);
    for (size_t i = 0, j = 0; i < volume_vox && j < target_nvox; ++i) {
        if (total_mask.at(i) > 0.5) {
            if (target_mask.data.at(i) > 0.5) {
                target_dose_vol[i] = target_dose_missing[j];
            }
            j++;
        }
    }
    Volume_t temp_vol1(target_mask.getMetadata());
    std::copy(std::begin(target_dose_vol), std::end(target_dose_vol), temp_vol1.data.begin());
    temp_vol1.output(out_directory+"/target_dose_missing.mhd");
    target_dose_vol.resize(0);
    Volume_t temp_vol2(target_mask.getMetadata());
    std::copy(std::begin(dose_by_subset), std::end(dose_by_subset), temp_vol2.data.begin());
    temp_vol2.output(out_directory+"/dose_in_subset_dij.mhd");
    dose_by_subset.resize(0);


    // 5: Prepare and launch Opt4D and get the weigth scaling ----------------
    Opt4D_manager opt4d(out_directory);
    opt4d.populate_directory(nspots_subset, dij_nvox, target_nvox, rim_nvox, oars_nvox,
                             total_mask, target_mask, target_rim_mask, oars_mask,
                             target_dose_missing, subset_dij);
    opt4d.launch_optimization();

    // 6: Apply weight scaling and set in referenced vector
    std::vector<float> subset_weight_scaling = opt4d.get_weight_scaling();

    std::vector<uint> accu_subset_per_field(nbeams);
    std::partial_sum(subset_per_field.begin(), subset_per_field.end(),
                     accu_subset_per_field.begin());

    for (size_t ibeam = 0; ibeam < indexes.size(); ++ibeam) {
        uint base_abs = 0;
        if (ibeam > 0)
            base_abs += accu_subset_per_field.at(ibeam-1);
        for (size_t i = 0; i < indexes.at(ibeam).size(); ++i) {
            // Set in referenced vector
            uint abs_spot = indexes.at(ibeam).at(i);
            if (ibeam > 0)
                abs_spot += accu_spots.at(ibeam-1);
            weight_scaling.at(abs_spot) = subset_weight_scaling.at(base_abs+i);
            
            // Scale dij subset dose
            Array temp = subset_dij.at(ibeam)[std::slice(i*dij_nvox, dij_nvox, 1)];
            temp *= subset_weight_scaling.at(base_abs+i);
            subset_dij.at(ibeam)[std::slice(i*dij_nvox, dij_nvox, 1)] = temp;
        }
    }

    // Add total dose by subset dij to baseline
    Array adapted_dose = baseline_dose;
    for (size_t i = 0, idij = 0; i < volume_vox && idij < dij_nvox; ++i) {
        if (total_mask.at(i) > 0.5) {
            float d = 0;
            for (size_t f = 0; f < nbeams; ++f) {
                Array temp = subset_dij[f][std::slice(idij, indexes.at(f).size(), dij_nvox)];
                d += temp.sum();
            }
            adapted_dose[i] += d;
            idij++;
        }
    }

    // 7: Check new plan dose
    check_adaptation(adapted_dose, dose_prescription,
                     target_nvox, oars_nvox, dij_nvox, total_mask, target_mask, oars_mask,
                     out_directory);
}


// void output_debug_doses (std::string out_directory, const Array& underdose_mask,
//                          const Volume_t& target_mask,
//                          const std::vector<std::string>& plan_dose_per_field_files,
//                          const Array& adapt_dose,
//                          const std::vector<Array>& adapt_field_dose,
//                          const Array& adapt_dose_in_target,
//                          float dose_prescription)
// {
//     std::ofstream ofs0(out_directory+"/adapt_dose.dat", std::ios::binary);
//     ofs0.write((char*)&adapt_dose[0], adapt_dose.size()*sizeof(float));
//     ofs0.close();

//     Array temp;
//     temp.resize(adapt_dose.size());
//     for (size_t i = 0, itarget = 0; i < adapt_dose.size(); ++i) {
//         if (target_mask.data.at(i) > 0.5) {
//             temp[i] = adapt_dose_in_target[itarget];
//             itarget++;
//         }
//     }
//     // std::ofstream ofs5(out_directory+"/adapt_dose_to_target.dat", std::ios::binary);
//     // ofs5.write((char*)&temp[0], temp.size()*sizeof(float));
//     // ofs5.close();

//     std::ofstream ofs(out_directory+"/cold_spots_map.dat", std::ios::binary);
//     ofs.write((char*)&underdose_mask[0], underdose_mask.size()*sizeof(float));
//     ofs.close();

//     // std::ofstream ofs2(out_directory+"/target_mask.dat", std::ios::binary);
//     // ofs2.write((char*)target_mask.data.data(), target_mask.data.size()*sizeof(float));
//     // ofs2.close();

//     // uint nbeams = plan_dose_per_field_files.size();
//     // uint volume_vox = target_mask.nElements;
//     // std::vector<std::vector<float>> diffs(nbeams, std::vector<float>(volume_vox, 0));
//     // Array total_plan_dose;
//     // total_plan_dose.resize(volume_vox);
//     // for (size_t i = 0; i < nbeams; ++i) {
//     //     Volume_t plan_dose(plan_dose_per_field_files.at(i), Volume_t::Source_type::DOSE);
//     //     for (size_t j = 0; j < volume_vox; ++j) {
//     //         total_plan_dose[j] += plan_dose.data.at(j);
//     //         if (target_mask.data.at(j) > 0.5 &&
//     //             adapt_dose[j] < dose_prescription &&
//     //             adapt_field_dose.at(i)[j] < plan_dose.data.at(j)) {
//     //             // Fill diff per beam with cold spot info
//     //             diffs.at(i).at(j) = plan_dose.data.at(j) - adapt_field_dose.at(i)[j];
//     //         }
//     //     }
//     //     std::ofstream ofs3(out_directory+"/cold_spots_map.beam_" +
//     //                        std::to_string(i)+".dat", std::ios::binary);
//     //     ofs3.write((char*)diffs.at(i).data(), diffs.at(i).size()*sizeof(float));
//     //     ofs3.close();
//     // }
// }
